      subroutine ringfilter(f,lev0,lev1,order,name,lat)

! Coded by Tong Dang, May, 2017; Modified by Haonan Wu, Feb, 2022
! This program is to build a new ring filter for the high-resolution TIE-GCM.
! The filter procedures are as follows:
! 1. Divide the longitude sessions into several chunks
! 2. Calculate the avg values of each chunk
! 3. Do the restruction based on the mean value of each chunk

      use params_module,only: nlon,nlonp2,nlonp4,nlat,rp
      use cons_module,only: dlamda,nlat_filter,chunk_array
      implicit none

      integer,intent(in) :: lev0,lev1,order,lat
      real(rp),dimension(nlonp4,lev0:lev1),intent(inout) :: f
      character(len=*),intent(in) :: name

      integer :: i,m,n,nchnk,ichnk,im2,im1,ip1,ip2
      real(rp) :: ang
      real(rp),dimension(lev0:lev1) :: a0,left,right,fa,fb
      real(rp),dimension(order,lev0:lev1) :: a,b
      real(rp),dimension(nlon,lev0:lev1) :: w,wm,fx,avg
      real(rp),dimension(nlon,lev0:lev1) :: spline_a, spline_b,
     |    spline_c, spline_d
      real(rp),dimension(nlon,lev0:lev1) :: poly_a,poly_b,poly_c,poly_d

      nchnk = 0
      if (lat <= nlat_filter) then
        nchnk = chunk_array(lat)
      elseif (lat >= nlat-nlat_filter+1) then
        nchnk = chunk_array(nlat-lat+1)
      else
        if (trim(name) == 'OPLUS') then
          i = max(abs(nlat/2-lat),1)
!          if (i <= nlat_filter) nchnk = chunk_array(i)
        endif
      endif
      if (nchnk == 0) return

      if (order >= nlon) call shutdown('Filtering order too high')

      w = f(3:nlonp2,:)

! Fourier series to nth order
      a0 = sum(w,1)/nlon
      a = 0._rp
      b = 0._rp
      do i = 1,nlon
        ang = i*dlamda
        do m = 1,order
          a(m,:) = a(m,:)+w(i,:)*cos(m*ang)
          b(m,:) = b(m,:)+w(i,:)*sin(m*ang)
        enddo
      enddo
      a = 2._rp*a/nlon
      b = 2._rp*b/nlon

! reconstruct based on Fourier series (taken as mean)
      do i = 1,nlon
        ang = i*dlamda
        wm(i,:) = a0
        do m = 1,order
          wm(i,:) = wm(i,:)+a(m,:)*cos(m*ang)+b(m,:)*sin(m*ang)
        enddo
      enddo

! deviation field
      fx = w-wm

! divide lon array into chunks and calculate mean in each chunk
      n = nlon/nchnk

      do ichnk = 1,nchnk
        avg(ichnk,:) = sum(fx((ichnk-1)*n+1:ichnk*n,:),1)/n
      enddo

      do ichnk = 1,nchnk
        im2 = ichnk-2
        if (im2 <= 0) im2 = im2+nchnk
        im1 = ichnk-1
        if (im1 <= 0) im1 = im1+nchnk
        ip1 = ichnk+1
        if (ip1 >= nchnk+1) ip1 = ip1-nchnk
        ip2 = ichnk+2
        if (ip2 >= nchnk+1) ip2 = ip2-nchnk
        left = (-avg(im2,:)+7._rp*avg(im1,:)+
     |  7._rp*avg(ichnk,:)-avg(ip1,:))/12._rp
        right = (-avg(im1,:)+7._rp*avg(ichnk,:)+
     |  7._rp*avg(ip1,:)-avg(ip2,:))/12._rp
        fa = (left+right-2._rp*avg(ichnk,:))/n**2
        fb = (3._rp*avg(ichnk,:)-2._rp*left-right)/n

        if (trim(name) == 'OPLUS') then
! piece-wise polynomial reconstruction
          do m = 1,n
            fx((ichnk-1)*n+m,:) = fa*(3._rp*m**2-3._rp*m+1._rp)+
     |  fb*(2._rp*m-1._rp)+left
          enddo
!
!            spline_a(ichnk,:) = avg(ichnk,:)
!            spline_b(ichnk,:) = (avg(ip1,:) - avg(im1,:)) / (2.0_rp * n)
!            spline_c(ichnk,:) = (avg(im1,:) - 2.0_rp*avg(ichnk,:) + avg
!     |  (ip1,:)) / (2.0_rp * n**2)
!            spline_d(ichnk,:) = (avg(ip1,:) - avg(ichnk,:)) / (n**3) -
!     |  spline_b(ichnk,:) / n
!
!            do m = 1,n
!            fx((ichnk-1)*n + m,:) = spline_a(ichnk,:) +
!     |                               spline_b(ichnk,:) * m +
!     |                               spline_c(ichnk,:) * m**2 +
!     |                               spline_d(ichnk,:) * m**3
!            enddo
!
!        ! Polynomial coefficients
!        poly_a(ichnk,:) = avg(ichnk,:)
!        poly_b(ichnk,:) = (avg(ip1,:) - avg(im1,:)) / (2.0_rp * n)
!        poly_c(ichnk,:) = (avg(im1,:) - 2.0_rp*avg(ichnk,:) + avg(ip1,:)
!     |     ) / (2.0_rp * n**2)
!        poly_d(ichnk,:) = (avg(ip1,:) - avg(ichnk,:)) / (n**3) - poly_b(
!     |     ichnk,:) / n
!        do m = 1,n
!            fx((ichnk-1)*n + m,:) = poly_a(ichnk,:) +
!     &                               poly_b(ichnk,:) * m +
!     &                               poly_c(ichnk,:) * m**2 +
!     &                               poly_d(ichnk,:) * m**3
!        enddo
!        left = (-avg(im2, :) + 7._rp * avg(im1, :) + 7._rp *
!     |    avg(ichnk, :) - avg(ip1, :)) / 12._rp
!        right = (-avg(im1, :) + 7._rp * avg(ichnk, :) + 7._rp *
!     |    avg(ip1, :) - avg(ip2, :)) / 12._rp
!        poly_a(ichnk, :) = avg(ichnk, :)
!        poly_b(ichnk, :) = (right - left) / (2.0_rp * n)
!        poly_c(ichnk, :) = (left - 2.0_rp * avg(ichnk, :) +
!     |    right) / (2.0_rp * n**2)
!        poly_d(ichnk, :) = (avg(ip1, :) - avg(ichnk, :)) /
!     |    (n**3) - poly_b(ichnk, :) / n
!        do m = 1, n
!            fx((ichnk-1)*n + m, :) = poly_a(ichnk, :) +
!     |           poly_b(ichnk, :) * m +
!     |           poly_c(ichnk, :) * m**2 +
!     |           poly_d(ichnk, :) * m**3
!        enddo
        else
! piece-wise linear reconstruction
          do m = 1,n
            fx((ichnk-1)*n+m,:) = (left*(n-m)+right*m)/n
          enddo
        endif
      enddo

! full field
      w = fx+wm

      f(3:nlonp2,:) = w
      f(1:2,:) = w(nlon-1:nlon,:)
      f(nlonp4-1:nlonp4,:) = w(1:2,:)

      end subroutine ringfilter
