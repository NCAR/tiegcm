module subaur_module
!-----------------------------------------------------------------------
! Copyright 2009 The Johns Hopkins University Applied Physics Laboratory
! All rights reserved.
!
! +
!   SAPS Electric Field Model
!   Principal Investigator and Contact:
!   Dr. Elsayed Talaat, Johns Hopkins University Applied Physics Laboratory
! +
!
! Software History:
!   Aug 27, 2009 - version 1.0
!   Sep 2020, refactored using F90 grammar, by Haonan Wu
!
! Author:
!   Syau-Yun Hsieh, Johns Hopkins University, Applied Physics Laboratory
!-----------------------------------------------------------------------
  implicit none

  integer, parameter :: nmlt=48, nalat=40, nkp=9
  real, dimension(nmlt) :: mlt
  real, dimension(nalat) :: alat
  real, dimension(nkp) :: kp
  real, dimension(nmlt, nalat, nkp) :: ve, vn, vv

  contains
!-----------------------------------------------------------------------
! Name:
!   init_saps
!
! Purpose:
!   - used to initialize the SAPS Model 
!
! Calling Sequence:
!   call init_saps
!
! Inputs:
!   None
!
! Outputs:
!   None
!-----------------------------------------------------------------------
  subroutine init_saps

    use input_module, only: vedata, vndata, vvdata

    integer :: ios, i, j, k
    external :: shutdown

    do i = 1, nkp
      kp(i) = i
    enddo
    do i = 1, nmlt
      mlt(i) = (i-1)*0.5
    enddo
    do i = 1, nalat
      alat(i) = (i-1)*0.25
    enddo

    open(unit=79, file=trim(vedata),iostat=ios,status='old')
    if (ios /= 0) call shutdown('Error opening VEDATA')
    read(79,*)(((ve(i,j,k),i=1,nmlt),j=1,nalat),k=1,nkp)
    close(79)

    open(unit=79, file=trim(vndata),iostat=ios,status='old')
    if (ios /= 0) call shutdown('Error opening VNDATA')
    read(79,*)(((vn(i,j,k),i=1,nmlt),j=1,nalat),k=1,nkp)
    close(79)

    open(unit=79, file=trim(vvdata),iostat=ios,status='old')
    if (ios /= 0) call shutdown('Error opening VVDATA')
    read(79,*)(((vv(i,j,k),i=1,nmlt),j=1,nalat),k=1,nkp)
    close(79)

  end subroutine init_saps
!-----------------------------------------------------------------------
! Name:
!   subaur_drift
!
! Purpose:
!   - used to obtain the subauroral drift [v-north,v-east,v-vertical]
!     from SAPS electric field model
!
! Calling Sequence:
!   call subaur_drift(ualat, umlt, ukp, vn1, ve1, vv1)
!
! Inputs:
!   ualat - auroral latitudes
!   umlt  - magnetic local time
!   ukp   - Kp index
!
! Outputs:
!   None
!-----------------------------------------------------------------------
  subroutine subaur_drift(ualat, umlt, ukp, vn1, ve1, vv1)

    real, intent(in) :: ualat, umlt, ukp
    real, intent(out) :: ve1, vn1, vv1

    integer :: k1, k2, k3
    real :: yp1, ypn, &
      ye1_new, yn1_new, yv1_new, &
      ye2_new, yn2_new, yv2_new, &
      ye3_new, yn3_new, yv3_new
    real, dimension(nkp) :: x1, y12, ye1, yn1, yv1
    real, dimension(nmlt) :: x2, y22, ye2, yn2, yv2
    real, dimension(nalat) :: x3, y32, ye3, yn3, yv3

    yp1 = 1.e30
    ypn = 1.e30

    do k1 = 1, nkp
      do k2 = 1, nmlt
        do k3 = 1, nalat
          x3(k3)  = alat(k3)
          ye3(k3) = ve(k2,k3,k1)
          yn3(k3) = vn(k2,k3,k1)
          yv3(k3) = vv(k2,k3,k1)
        enddo

        call spline(x3,ye3,nalat,yp1,ypn,y32)
        call splint(x3,ye3,y32,nalat,ualat,ye3_new)
        x2(k2) = mlt(k2)
        ye2(k2) = ye3_new

        call spline(x3,yn3,nalat,yp1,ypn,y32)
        call splint(x3,yn3,y32,nalat,ualat,yn3_new)
        yn2(k2) = yn3_new

        call spline(x3,yv3,nalat,yp1,ypn,y32)
        call splint(x3,yv3,y32,nalat,ualat,yv3_new)
        yv2(k2) = yv3_new

      enddo

      call spline(x2,ye2,nmlt,yp1,ypn,y22) 
      call splint(x2,ye2,y22,nmlt,umlt,ye2_new)
      x1(k1) = kp(k1)
      ye1(k1) = ye2_new

      call spline(x2,yn2,nmlt,yp1,ypn,y22)
      call splint(x2,yn2,y22,nmlt,umlt,yn2_new)
      yn1(k1) = yn2_new

      call spline(x2,yv2,nmlt,yp1,ypn,y22)
      call splint(x2,yv2,y22,nmlt,umlt,yv2_new)
      yv1(k1) = yv2_new

    enddo

    call spline(x1,ye1,nkp,yp1,ypn,y12)
    call splint(x1,ye1,y12,nkp,ukp,ye1_new)

    call spline(x1,yn1,nkp,yp1,ypn,y12)
    call splint(x1,yn1,y12,nkp,ukp,yn1_new)

    call spline(x1,yv1,nkp,yp1,ypn,y12)
    call splint(x1,yv1,y12,nkp,ukp,yv1_new)

    ve1 = ye1_new
    vn1 = yn1_new
    vv1 = yv1_new

  end subroutine subaur_drift
!----------------------------------------------------------------------------
!   subroutine spline
!   - used for cubic spline interpolation
!   - adopted from Numerical Recipes
!----------------------------------------------------------------------------
  subroutine spline(x, y, n, yp1, ypn, y2)

    integer, intent(in) :: n
    real, intent(in) :: yp1, ypn
    real, intent(in), dimension(n) :: x, y
    real, intent(out), dimension(n) :: y2

    integer :: i
    real :: p, qn, sig, un
    real, dimension(500) :: u

    if (yp1 > 0.99e30) then
      y2(1) = 0.
      u(1)  = 0.
    else
      y2(1) = -0.5
      u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2,n-1
      sig   = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
      p     = sig * y2(i-1) + 2.
      y2(i) = (sig-1.) / p
      u(i)  = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo

    if (ypn > 0.99e30) then
      qn = 0.
      un = 0.
    else
      qn = 0.5
      un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do i = n-1,1,-1
      y2(i) = y2(i)*y2(i+1)+u(i)
    enddo

  end subroutine spline
!-----------------------------------------------------------------------
!   subroutine splint
!   - used to call subroutine spline for interpolations
!   - adopted from Numerical Recipes
!-----------------------------------------------------------------------
  subroutine splint(xa, ya, y2a, n, x, y)

    integer, intent(in) :: n
    real, intent(in) :: x
    real, intent(in), dimension(n) :: xa, ya, y2a
    real, intent(out) :: y

    integer :: k, khi, klo
    real :: a, b, h

    klo = 1
    khi = n
    do while (khi-klo > 1)
      k = (khi+klo)/2
      if (xa(k) > x) then
        khi = k
      else
        klo = k
      endif
    enddo

    h = xa(khi) - xa(klo)
!   if (h == 0.) write(6,*) 'bad xa input in splint'
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+ &
      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

  end subroutine splint
!-----------------------------------------------------------------------
end module subaur_module
