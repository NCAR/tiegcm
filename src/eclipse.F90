module eclipse_module

  implicit none

  integer,parameter :: nelement = 28, ncircum = 40, nevent = 155, nc2 = 315, nc3 = 340, maxiter = 50
  real,parameter :: PI = 3.14151926, dtr = PI/180, rtd = 180/PI, eps = 1e-6

  logical :: has_eclipse

! Observer constants
!   (0) North Latitude (degree)
!   (1) East Longitude (degree)
!   (2) Altitude (meter)
!   (3) West time zone (hour)
!   (4) rho sin O'
!   (5) rho cos O'
!   (6) index into the elements array for the eclipse in question (not used)
  real,dimension(0:5) :: observer_constants

! Eclipse Elements
!   (0) Julian date of maximum eclipse
!   (1) t0 - the TDT hour at which t=0
!   (2) tmin - the lowest allowed value of t
!   (3) tmax - the highest allowed value of t
!   (4) dUTC - the difference between the civilian "GMT" timescale and TDT
!   (5) dT - the difference between ut and TDT
!   (6-9) X0, X1, X2, X3 - X elements
!   (10-13) Y0, Y1, Y2, Y3 - Y elements
!   (14-16) D0, D1, D2 - D elements
!   (17-19) M0, M1, M2 - mu elements
!   (20-22) L10, L11, L12 - L1 elements
!   (23-25) L20, L21, L22 - L2 elements
!   (26) tan f1
!   (27) tan f2
  real,dimension(0:nelement-1) :: elements

! Eclipse circumstances
!   (0) Event type (C1=-2, C2=-1, Mid=0, C3=1, C4=2)
!   (1) t
! -- time-only dependent circumstances (and their per-hour derivatives) follow --
!   (2) x  (3) y  (4)  d (5) sin d  (6) cos d  (7) mu  (8) l1  (9) l2  (10) dx
!   (11) dy (12) dd  (13) dmu  (14) dl1  (15) dl2
! -- time and location dependent circumstances follow --
!   (16) h  (17) sin h  (18) cos h  (19) xi  (20) eta  (21) zeta  (22) dxi  (23) deta  (24) u
!   (25) v  (26) a  (27) b  (28) l1' (29) l2' (30) n^2
! -- observational circumstances follow --
!   (31) p  (32) alt  (33) q  (34) v  (35) azi
!   (36) m (mid eclipse only) or limb correction applied (where available!)
!   (37) magnitude (mid eclipse only)
!   (38) moon/sun (mid eclipse only)
!   (39) calculated local event type for a transparent earth (0 = none, 1 = partial, 2 = annular, 3 = total)
  real,dimension(0:ncircum-1) :: c1,c2,c3,c4,mid

  integer,dimension(nevent,3) :: date_list
  real,dimension(nevent,nelement) :: eclipse_list
  real,dimension(0:nc2-1) :: C2limb2003May
  real,dimension(0:nc3-1) :: C3limb2003May

  contains
!-----------------------------------------------------------------------
  subroutine read_eclipse_list(filename)
! Read eclipse list from external data, called once per run

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_nowrite,nf90_noerr

    character(len=*),intent(in) :: filename

    integer :: stat,ncid,dimid,varid,n

    stat = nf90_open(trim(filename),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_open',stat)

    stat = nf90_inq_dimid(ncid,'event',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=n)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    if (n /= nevent) write(6,"('Eclipse data file dimension event does not conform')")

    stat = nf90_inq_dimid(ncid,'element',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=n)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    if (n /= nelement) write(6,"('Eclipse data file dimension element does not conform')")

    stat = nf90_inq_dimid(ncid,'c2',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=n)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    if (n /= nc2) write(6,"('Eclipse data file dimension c2 does not conform')")

    stat = nf90_inq_dimid(ncid,'c3',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=n)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    if (n /= nc3) write(6,"('Eclipse data file dimension c3 does not conform')")

    stat = nf90_inq_varid(ncid,'date',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,date_list)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_inq_varid(ncid,'eclipse',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,eclipse_list)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_inq_varid(ncid,'C2limb2003May',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,C2limb2003May)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_inq_varid(ncid,'C3limb2003May',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,C3limb2003May)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

  end subroutine read_eclipse_list
!-----------------------------------------------------------------------
  subroutine find_eclipse(year,month,day)
! Set the elements array

    integer,intent(in) :: year,month,day

    integer :: ievent

    do ievent = 1,nevent
      if (date_list(ievent,1)==year .and. &
        date_list(ievent,2)==month .and. &
        date_list(ievent,3)==day) exit
    enddo

    if (ievent <= nevent) then
      has_eclipse = .true.
      elements = eclipse_list(ievent,:)
    else ! There is no eclipse at this time
      has_eclipse = .false.
    endif

  end subroutine find_eclipse
!-----------------------------------------------------------------------
  subroutine eclipse(year,month,day,latitude,longitude,height,ut,slt,mag,fac)
! This subroutine must be called after find_eclipse

    integer,intent(in) :: year,month,day
    real,intent(in) :: latitude,longitude,height,ut,slt
    real,intent(out) :: mag,fac

    real :: t1,t4,tm,tshift,t1_slt,t4_slt,R,M,D,d1,S1, &
! For calculating the ratio of phi_Ti/phi_T0 = (12.3*Sp+Sc)/(12.3*Sp0+Sc0)
      Rp,Rc,Sp0,Sc0,Sp,Sc,dl,E,Xd,A1,A2,A3

    if (has_eclipse) then
      call set_observer_constants(latitude,longitude,height)
      call set_circumstances(year,month,day)

      call get_time(c1,t1)
      call get_time(c4,t4)
      call get_time(mid,tm)

      tshift = observer_constants(1)/15
      t1_slt = modulo(t1+tshift, 24.0)
      t4_slt = modulo(t4+tshift, 24.0)

! Calculate the actual ratio of the obscured diameter of the solar photosphere
      R = 0.5
      M = floor(100000*mid(37)+0.5)/100000
      D = sqrt(4*M*R-M**2)
      if (ut <= tm) then 
        d1 = (tm-ut)/(tm-t1)*D
      else
        d1 = (tm-ut)/(tm-t4)*D
      endif
      S1 = sqrt(d1**2+(2*R-M)**2)
      mag = 2*R-S1

      if (mid(39)==0 .or. & ! there is no eclipse
        ((t1_slt< t4_slt .and. (slt<t1_slt .or.  slt>t4_slt)) .or. &
         (t1_slt>=t4_slt .and. (slt>t4_slt .and. slt<t1_slt)))) then
        mag = 0
      endif

! Calculate the unmasked area for photosphere for P and C regions during the solar eclipse.
! P region is an inner central disk over the "photosphere";
! C region is an outer ring between the "coronal" and "phosphere".
! reference to J. J. Curto et al. (2006),
! Modeling the geomagnetic effects caused by the solar eclipse of 11 August 1999

! Initialize
      Rp = 1
      Rc = 2.1
      Sp0 = PI*Rp**2
      Sc0 = PI*(Rc**2-Rp**2)

      dl = 2-mag*2
      if (dl < Rc-Rp) then 
        A1 = 0
        E = 2*Rp**2*acos(dl/2/Rp) - dl*sqrt(Rp**2 - dl**2/4)
      elseif (dl < 2*Rp) then
        E = 2*Rp**2*acos(dl/2/Rp) - dl*sqrt(Rp**2 - dl**2/4)
        Xd = (Rc**2 - Rp**2 + dl**2)/2/dl
        A2 = Rc**2*acos(Xd/Rc) - Xd*sqrt(Rc**2 - Xd**2)
        A3 = Rp**2*acos((dl-Xd)/Rp) - (dl-Xd)*sqrt(Rp**2-(dl-Xd)**2)
        A1 = PI*Rp**2 - A2 - A3
      elseif (dl < Rc+Rp) then
        E = 0
        Xd = (Rc**2 - Rp**2 + dl**2)/2/dl
        A2 = Rc**2*acos(Xd/Rc) - Xd*sqrt(Rc**2 - Xd**2)
        A3 = Rp**2*acos((dl-Xd)/Rp) - (dl-Xd)*sqrt(Rp**2-(dl-Xd)**2)
        A1 = PI*Rp**2 - A2 - A3
      else
        write(6,"('dl should be smaller than Rp+Rc=3.1 while it is ',e12.5)") dl
      endif

      Sc = PI*(Rc**2 - 2*Rp**2) + E + A1
      Sp = PI*Rp**2 - E

      if (mag < 1e-4) then
        fac = 1
      else
        fac = (12.3*Sp+Sc)/(12.3*Sp0+Sc0)
      endif
    else ! No eclipse
      mag = 0
      fac = 1
    endif

  end subroutine eclipse
!-----------------------------------------------------------------------
  subroutine set_observer_constants(latitude,longitude,height)

    real,intent(in) :: latitude,longitude,height

    real,parameter :: aspect_ratio = 0.99664719, re = 6378140
    real :: latr,reduced_lat

! Note that correcting for refraction will involve creating a "virtual" altitude
! for each contact, and hence a different value of rho and O' for each contact!
! There was a minus sign in front of longitude, change to plus throughout the module
    observer_constants(0) = latitude
    observer_constants(1) = longitude
    observer_constants(2) = height*1e3
    observer_constants(3) = 0

    latr = observer_constants(0)*dtr
    reduced_lat = atan(aspect_ratio*tan(latr))
    observer_constants(4) = aspect_ratio*sin(reduced_lat) + observer_constants(2)/re*sin(latr)
    observer_constants(5) =              cos(reduced_lat) + observer_constants(2)/re*cos(latr)

  end subroutine set_observer_constants
!-----------------------------------------------------------------------
  subroutine set_circumstances(year,month,day)
! Populate the c1, c2, mid, c3 and c4 arrays

    integer,intent(in) :: year,month,day

    call set_mid
    call observational(mid)

! Calculate m, magnitude and moon/sun
    mid(36) = sqrt(mid(24)**2 + mid(25)**2)
    mid(37) = (mid(28)-mid(36)) / (mid(28)+mid(29))
    mid(38) = (mid(28)-mid(29)) / (mid(28)+mid(29))
    if (mid(37) > 0) then
      call set_c1c4
      if (mid(36)<mid(29) .or. mid(36)<-mid(29)) then
        call set_c2c3
        if (mid(29) < 0) then
          mid(39) = 3 ! Total eclipse
        else 
          mid(39) = 2 ! Annular eclipse
        endif
        call observational(c2)
        call observational(c3)

! 2003 May 31 eclipse limb corrections
        if (year==2003 .and. month==5 .and. day==31) then
          call limb_correction(c2(31),c3(31),c2(36),c3(36))
          if (c2(36) < 990) c2(1) = c2(1) + c2(36)/3600
          if (c3(36) < 990) c3(1) = c3(1) + c3(36)/3600
        else 
          c2(36) = 999.9
          c3(36) = 999.9
        endif
      else 
        mid(39) = 1 ! Partial eclipse
      endif
      call observational(c1)
      call observational(c4)
    else 
      mid(39) = 0 ! No eclipse
    endif

  end subroutine set_circumstances
!-----------------------------------------------------------------------
  subroutine set_mid
! Calculate mid eclipse

    integer :: iter
    real :: tmp

    mid(0) = 0
    mid(1) = 0
    call calculate_circumstances(mid)

    tmp = 1
    do iter = 1,maxiter
      if (abs(tmp) <= eps) exit

      tmp = (mid(24)*mid(26) + mid(25)*mid(27)) / mid(30)
      mid(1) = mid(1) - tmp
      call calculate_circumstances(mid)
    enddo

  end subroutine set_mid
!-----------------------------------------------------------------------
  subroutine observational(circumstances)
! Get the observational circumstances

    real,dimension(0:ncircum-1),intent(inout) :: circumstances

    integer :: typepe
    real :: contacttype,coslat,sinlat

    typepe = circumstances(0)

! We are looking at an "external" contact UNLESS this is a total eclipse AND we are looking at
! c2 or c3, in which case it is an INTERNAL contact! Note that if we are looking at mid eclipse,
! then we may not have determined the type of eclipse (mid[39]) just yet!
    if (mid(39)==3 .and. abs(typepe)==1) then
      contacttype = -1
    else
      contacttype = 1
    endif

! Calculate p
    circumstances(31) = atan2(contacttype*circumstances(24), contacttype*circumstances(25))

! Calculate altitude
    sinlat = sin(observer_constants(0)*dtr)
    coslat = cos(observer_constants(0)*dtr)
    circumstances(32) = asin(circumstances(5)*sinlat + circumstances(6)*coslat*circumstances(18))

! Calculate q
    circumstances(33) = asin(coslat*circumstances(17)/cos(circumstances(32)))
    if (circumstances(20) < 0) circumstances(33) = PI-circumstances(33)

! Calculate v
    circumstances(34) = circumstances(31)-circumstances(33)

! Calculate azimuth
    circumstances(35) = atan2(-circumstances(17)*circumstances(6), &
      circumstances(5)*coslat-circumstances(18)*sinlat*circumstances(6))

  end subroutine observational
!-----------------------------------------------------------------------
  subroutine set_c1c4
! Set C1 and C4 data

! Entry conditions -
!   1. The mid array must be populated
!   2. The magnitude at mid eclipse must be > 0

    real :: tmp,n

    n = sqrt(mid(30))
    tmp = (mid(26)*mid(25) - mid(24)*mid(27)) / n / mid(28)
    tmp = sqrt(max(1-tmp**2,0.0))*mid(28)/n
    c1(0) = -2
    c4(0) = 2
    c1(1) = mid(1)-tmp
    c4(1) = mid(1)+tmp
    call iterate_c1c4(c1)
    call iterate_c1c4(c4)

  end subroutine set_c1c4
!-----------------------------------------------------------------------
  subroutine iterate_c1c4(circumstances)
! Iterate on C1 or C4

    real,dimension(0:ncircum-1),intent(inout) :: circumstances

    integer :: iter
    real :: sign,tmp,n

    call calculate_circumstances(circumstances)
    if (circumstances(0) < 0) then
      sign = -1
    else 
      sign = 1
    endif

    tmp = 1
    do iter = 1,maxiter
      if (abs(tmp) <= eps) exit

      n = sqrt(circumstances(30))
      tmp = (circumstances(26)*circumstances(25) - circumstances(24)*circumstances(27)) / n / circumstances(28)
      tmp = sign*sqrt(max(1-tmp**2,0.0))*circumstances(28)/n
      tmp = (circumstances(24)*circumstances(26) + circumstances(25)*circumstances(27)) / circumstances(30) - tmp
      circumstances(1) = circumstances(1)-tmp
      call calculate_circumstances(circumstances)
    enddo

  end subroutine iterate_c1c4
!-----------------------------------------------------------------------
  subroutine set_c2c3
! Set C2 and C3 data

! Entry conditions -
!   1. The mid array must be populated
!   2. There mut be either a total or annular eclipse at the location!

    real :: tmp,n
  
    n = sqrt(mid(30))
    tmp = (mid(26)*mid(25) - mid(24)*mid(27)) / n / mid(29)
    tmp = sqrt(max(1-tmp**2,0.0))*mid(29)/n
    c2(0) = -1
    c3(0) = 1
    if (mid(29) < 0) then
      c2(1) = mid(1)+tmp
      c3(1) = mid(1)-tmp
    else 
      c2(1) = mid(1)-tmp
      c3(1) = mid(1)+tmp
    endif
    call iterate_c2c3(c2)
    call iterate_c2c3(c3)

  end subroutine set_c2c3
!-----------------------------------------------------------------------
  subroutine iterate_c2c3(circumstances)
! Iterate on C2 or C3

    real,dimension(0:ncircum-1),intent(inout) :: circumstances

    integer :: iter
    real :: sign,tmp,n

    call calculate_circumstances(circumstances)
    if (circumstances(0) < 0) then
      sign = -1
    else 
      sign = 1
    endif
    if (mid(29) < 0) sign = -sign

    tmp = 1
    do iter = 1,maxiter
      if (abs(tmp) <= eps) exit

      n = sqrt(circumstances(30))
      tmp = (circumstances(26)*circumstances(25) - circumstances(24)*circumstances(27)) / n / circumstances(29)
      tmp = sign*sqrt(max(1-tmp**2,0.0))*circumstances(29)/n
      tmp = (circumstances(24)*circumstances(26) + circumstances(25)*circumstances(27)) / circumstances(30) - tmp
      circumstances(1) = circumstances(1)-tmp
      call calculate_circumstances(circumstances)
    end do

  end subroutine iterate_c2c3
!-----------------------------------------------------------------------
  subroutine calculate_circumstances(circumstances)
! Populate the circumstances array with the time and location dependent circumstances

    real,dimension(0:ncircum-1),intent(inout) :: circumstances

    integer :: typepe
    real :: t

    typepe = circumstances(0)
    t = circumstances(1)

! Calculate x
    circumstances(2) = ((elements(9)*t+elements(8))*t+elements(7))*t+elements(6)

! Calculate y
    circumstances(3) = ((elements(13)*t+elements(12))*t+elements(11))*t+elements(10)

! Calculate d
    circumstances(4) = ((elements(16)*t+elements(15))*t+elements(14))*dtr

! sin d and cos d
    circumstances(5) = sin(circumstances(4))
    circumstances(6) = cos(circumstances(4))

! Calculate m
    circumstances(7) = modulo((elements(19)*t+elements(18))*t+elements(17), 360.0) * dtr

! Calculate l1, dl1 and l1'
    if (typepe==-2 .or. typepe==0 .or. typepe==2) then
      circumstances(8) = (elements(22)*t+elements(21))*t+elements(20)
      circumstances(14) = 2*elements(22)*t+elements(21)
      circumstances(28) = circumstances(8)-circumstances(21)*elements(26)
    endif

! Calculate l2, dl2 and l2'
    if (typepe==-1 .or. typepe==0 .or. typepe==1) then
      circumstances(9) = (elements(25)*t+elements(24))*t+elements(23)
      circumstances(15) = 2*elements(25)*t+elements(24)
      circumstances(29) = circumstances(9)-circumstances(21)*elements(27)
    endif

! Calculate dx
    circumstances(10) = (3*elements(9)*t+2*elements(8))*t+elements(7)

! Calculate dy
    circumstances(11) = (3*elements(13)*t+2*elements(12))*t+elements(11)

! Calculate dd
    circumstances(12) = (2*elements(16)*t+elements(15))*dtr

! Calculate dm
    circumstances(13) = (2*elements(19)*t+elements(18))*dtr

! Calculate h, sin h, cos h
    circumstances(16) = circumstances(7) + observer_constants(1)*dtr - elements(5)/13713.44
    circumstances(17) = sin(circumstances(16))
    circumstances(18) = cos(circumstances(16))

! Calculate xi
    circumstances(19) = observer_constants(5)*circumstances(17)

! Calculate eta
    circumstances(20) = observer_constants(4)*circumstances(6) - observer_constants(5)*circumstances(18)*circumstances(5)

! Calculate zeta
    circumstances(21) = observer_constants(4)*circumstances(5) + observer_constants(5)*circumstances(18)*circumstances(6)

! Calculate dxi
    circumstances(22) = circumstances(13)*observer_constants(5)*circumstances(18)

! Calculate deta
    circumstances(23) = circumstances(13)*circumstances(19)*circumstances(5) - circumstances(21)*circumstances(12)

! Calculate u
    circumstances(24) = circumstances(2)-circumstances(19)

! Calculate v
    circumstances(25) = circumstances(3)-circumstances(20)

! Calculate a
    circumstances(26) = circumstances(10)-circumstances(22)

! Calculate b
    circumstances(27) = circumstances(11)-circumstances(23)

! Calculate n^2
    circumstances(30) = circumstances(26)**2 + circumstances(27)**2

  end subroutine calculate_circumstances
!-----------------------------------------------------------------------
  subroutine get_time(circumstances,t)
! Get the local time of an event

    real,dimension(0:ncircum-1),intent(in) :: circumstances
    real,intent(out) :: t

! Calculate the local time. Add 0.05 seconds, as we will be rounding up to the nearest 0.1 sec
    t = modulo(circumstances(1)+elements(1)-observer_constants(3)-(elements(4)-0.05)/3600, 24.0)

  end subroutine get_time
!-----------------------------------------------------------------------
  subroutine get_altitude(circumstances,alt)

    real,dimension(0:ncircum-1),intent(in) :: circumstances
    real,intent(out) :: alt

    alt = circumstances(32)*rtd

  end subroutine get_altitude
!-----------------------------------------------------------------------
  subroutine get_azimuth(circumstances,azi)

    real,dimension(0:ncircum-1),intent(in) :: circumstances
    real,intent(out) :: azi

    azi = modulo(circumstances(35)*rtd, 360.0)

  end subroutine get_azimuth
!-----------------------------------------------------------------------
  subroutine get_p(circumstances,p)

    real,dimension(0:ncircum-1),intent(in) :: circumstances
    real,intent(out) :: p

    p = modulo(circumstances(31)*rtd, 360.0)

  end subroutine get_p
!-----------------------------------------------------------------------
  subroutine get_v(circumstances,v)

    real,dimension(0:ncircum-1),intent(in) :: circumstances
    real,intent(out) :: v

    v = floor(120.5 - circumstances(34)*60/PI) / 10
    do while (v > 13)
      v = v-12
    enddo
    if (v < 1) v = v+12

  end subroutine get_v
!-----------------------------------------------------------------------
  subroutine get_duration(duration)

    real,intent(out) :: duration

    duration = modulo(c3(1)-c2(1), 24.0)

  end subroutine get_duration
!-----------------------------------------------------------------------
  subroutine get_coverage(cover)

    real,intent(out) :: cover

    real :: a,b,c
  
    if (mid(37) <= 0) cover = 0
    if (mid(37) >= 1) cover = 1
    if (mid(39) == 2) then
      c = mid(38)**2
    else
      c = acos((mid(28)**2 + mid(29)**2 - 2*mid(36)**2) / (mid(28)**2 - mid(29)**2))
      b = acos((mid(28)*mid(29) + mid(36)**2) / mid(36) / (mid(28)+mid(29)))
      a = PI-b-c
      c = (mid(38)**2*a + b - mid(38)*sin(c))/PI
    end if
    cover = c+0.5/1000

  end subroutine get_coverage
!-----------------------------------------------------------------------
  subroutine limb_correction(p2,p3,q2,q3)
! C2 and C3 limb corrections for the 2003 May 31 annular eclipse in seconds

! The data starts at contact angle 197.37/353.37 degrees in 0.4 degree increments for
! 315/340 data points, which means that the last element is for angle 322.97/128.97 degrees

! These limb corrections were calculated by Fred Espenak, NASA/GSFC

    real,intent(in) :: p2,p3
    real,intent(out) :: q2,q3

    real,parameter :: C2limb_start = 197.37, C3limb_start = 353.37, increment = 0.4
    integer :: n
    real :: frac

    frac = modulo(p2*rtd-C2limb_start, 360.0) / increment
    n = frac

    if (n >= nc2) then
      q2 = 999
    else
      q2 = (C2limb2003May(n+1)-C2limb2003May(n))*(frac-n) + C2limb2003May(n)
    endif

    frac = modulo(p3*rtd-C3limb_start, 360.0) / increment
    n = frac

    if (n >= nc3) then
      q3 = 999
    else
      q3 = (C3limb2003May(n+1)-C3limb2003May(n))*(frac-n) + C3limb2003May(n)
    endif

  end subroutine limb_correction
!-----------------------------------------------------------------------
  subroutine handle_error(funcname, ncerr)

    use netcdf, only: nf90_strerror

    character(len=*), intent(in) :: funcname
    integer, intent(in) :: ncerr

    write(6, "('NetCDF error encountered: ', a, ', when calling ', a)") &
      trim(nf90_strerror(ncerr)), funcname

  endsubroutine handle_error
!-----------------------------------------------------------------------
end module eclipse_module
