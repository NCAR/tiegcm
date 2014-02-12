! author Douglas P. Drob, Space Science Division, Naval Research Laboratory

module oxfmod

  implicit none

  real(8),parameter    :: deg2rad = 1.74533E-2
  real(8),parameter    :: day2rad = 1.72142E-2
  
  integer,parameter    :: nbf = 24
  integer,parameter    :: nalt = 19 
  real(8)              :: oxfp(0:nbf-1,0:nalt-1)
  real(8),parameter    :: altmin = 82.0
  real(8),parameter    :: dalt = 1.0
  logical              :: initparms = .true.
  
  
contains

  subroutine getparms()
    implicit none
    open(unit=1,file='oxparms.bin',access='stream')
    read(1) oxfp
    close(1)
    initparms = .false.
    return
  end subroutine getparms

  subroutine oxf(day,lat,alt,o)

    real(8),intent(in)   :: day
    real(8),intent(in)   :: lat
    real(8),intent(in)   :: alt
    real(8),intent(out)  :: o

    real(8)              :: slat,clat,cday,sday,c2day,s2day
    real(8)              :: legp2,legp3,legp4,legp5,legp6,legp7,legp8
    real(8)              :: f(0:nbf-1)
    real(8)              :: a,b
    integer              :: iz
    
    if (initparms) call getparms()

    slat = cos(deg2rad*lat)  ! slat <=> cos, legendre polynomial defined in colat
    clat = sin(deg2rad*lat)
    cday = cos(day2rad*day)
    sday = sin(day2rad*day)
    c2day = cos(day2rad*day*2)
    s2day = sin(day2rad*day*2)

    legp2 = clat
    legp3 = 0.5*(3.0*clat**2 - 1.0)
    legp4 = 0.5*(5.0*clat**3 - 3.0*clat)
    legp5 = (35.0*clat**4 - 30.0*clat**2 + 3.0)/8.0
    legp6 = (63.0*clat**5 - 70.*clat**3 + 15.0*clat)/8.0
    legp7 = (11.0*clat*legp6 - 5.0*legp5)/6.0
    legp8 = (13.0*clat*legp7 - 6.0*legp6)/7.0

    
    f(0) = 1.
    f(1) = cday
    f(2) = sday
    f(3) = cday*legp2
    f(4) = sday*legp2
    f(5) = cday*legp4
    f(6) = sday*legp4  
    f(7) = c2day
    f(8) = s2day
    f(9) = c2day*legp3
    f(10) = s2day*legp3
    f(11) = c2day*legp2
    f(12) = s2day*legp2  
    f(13) = legp3
    f(14) = legp5
    f(15) = legp7
    f(16) = legp2
    f(17) = cday*legp6
    f(18) = sday*legp6
    f(19) = c2day*legp5
    f(20) = s2day*legp5
    f(21) = cday*legp3
    f(22) = sday*legp3
    f(23) = legp4

    iz = int((alt - altmin)/dalt)
    if (iz .lt. 0) iz = 0
    if (iz .gt. nalt-2) iz = nalt-2
    b = (alt - (altmin + float(iz)*dalt))/dalt 
    a = 1.0 - b
    o = a*dot_product(oxfp(:,iz),f) + b*dot_product(oxfp(:,iz+1),f)
    
    return

  end subroutine oxf

end module oxfmod
