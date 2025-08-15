subroutine divergence(f,df,dir,lev0,lev1,lon0,lon1,lat0,lat1)

  implicit none

  integer,intent(in) :: dir,lev0,lev1,lon0,lon1,lat0,lat1
  real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in) :: f
  real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: df

  integer :: k,i,lat

  do lat = lat0,lat1
    do i = lon0,lon1
      do k = lev0,lev1
        if (dir == 1) then
          df(k,i,lat) = &
            2*(f(k,i,lat+1)-f(k,i,lat-1))/3- &
              (f(k,i,lat+2)-f(k,i,lat-2))/12
        endif
        if (dir == -1) then
          df(k,i,lat) = &
            2*(f(k,i+1,lat)-f(k,i-1,lat))/3- &
              (f(k,i+2,lat)-f(k,i-2,lat))/12
        endif
      enddo
    enddo
  enddo

endsubroutine divergence
