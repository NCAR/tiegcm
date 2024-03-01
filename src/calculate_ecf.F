      subroutine calculate_ecf(mag_factor,eclipse_factor,
     |  lon0,lon1,lev0,lev1,lat)
! calculate the rest EUV flux during eclipse period, Dang 2017

      use params_module,only: nlon,nlonp4,glat,glon
      use cons_module,only: re,pi
      use hist_module,only: modeltime
      use chapman_module,only: chi
      use init_module,only: iyear,iday
      use eclipse_module,only: eclipse

      implicit none

      integer,intent(in) :: lon0,lon1,lev0,lev1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(out) ::
     |  mag_factor,eclipse_factor

      integer :: month,day,i,ishift
      real :: height,rekm,sza,longitude,latitude,ut,slt,mag,fac

! Currently height is fixed at 300km for all altitudes
      height = 300

      rekm = re*1e-5
      sza = pi-asin(rekm/(rekm+height))

      call to_month_day(iyear,iday,month,day)
      ut = modeltime(2)+modeltime(3)/60+modeltime(4)/3600

      latitude = glat(lat)
      do i = lon0,lon1
        if (i==1 .or. i==2) then
          ishift = i-2+nlon
        elseif (i==nlonp4-1 .or. i==nlonp4) then
          ishift = i-2-nlon
        else
          ishift = i-2
        endif
        longitude = glon(ishift)
        slt = modulo(ut+longitude/15, 24.0)
        if (chi(i,lat) <= sza) then
          call eclipse(iyear,month,day,
     |      latitude,longitude,height,ut,slt,mag,fac)
          mag_factor(lev0:lev1,i) = mag
          eclipse_factor(lev0:lev1,i) = fac
        else
          mag_factor(lev0:lev1,i) = 0
          eclipse_factor(lev0:lev1,i) = 1
        endif
      enddo

      end subroutine calculate_ecf
