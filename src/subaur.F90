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

  integer :: nmlt, nalat, nkp
  real, dimension(:), allocatable :: mlt, alat, kp
  real, dimension(:, :, :), allocatable :: vn, ve, vv

  contains
!-----------------------------------------------------------------------
  subroutine init_saps(subaur_data)
! Name:
!   init_saps
!
! Purpose:
!   - used to initialize the SAPS Model 
!
! Calling Sequence:
!   call init_saps(subaur_data)
!
! Inputs:
!   subaur_data: filename of the SAPS data file
!
! Outputs:
!   None

    use netcdf, only: nf90_open, nf90_inq_dimid, nf90_inquire_dimension, &
      nf90_inq_varid, nf90_get_var, nf90_close, nf90_nowrite, nf90_noerr

    character(len=*), intent(in) :: subaur_data

    integer :: stat, ncid, dimid, varid

    stat = nf90_open(subaur_data, nf90_nowrite, ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_open', stat)

    stat = nf90_inq_dimid(ncid, 'mlt', dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid', stat)

    stat = nf90_inquire_dimension(ncid, dimid, len=nmlt)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension', stat)

    stat = nf90_inq_dimid(ncid, 'alat', dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid', stat)

    stat = nf90_inquire_dimension(ncid, dimid, len=nalat)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension', stat)

    stat = nf90_inq_dimid(ncid, 'kp', dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid', stat)

    stat = nf90_inquire_dimension(ncid, dimid, len=nkp)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension', stat)

    allocate(mlt(nmlt))
    allocate(alat(nalat))
    allocate(kp(nkp))
    allocate(vn(nmlt, nalat, nkp))
    allocate(ve(nmlt, nalat, nkp))
    allocate(vv(nmlt, nalat, nkp))

    stat = nf90_inq_varid(ncid, 'mlt', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, mlt)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_inq_varid(ncid, 'alat', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, alat)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_inq_varid(ncid, 'kp', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, kp)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_inq_varid(ncid, 'vn', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, vn)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_inq_varid(ncid, 've', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, ve)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_inq_varid(ncid, 'vv', varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid', stat)

    stat = nf90_get_var(ncid, varid, vv)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var', stat)

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close', stat)

  end subroutine init_saps
!-----------------------------------------------------------------------
  subroutine subaur_drift(umlt, ualat, ukp, vn1, ve1, vv1)
! Name:
!   subaur_drift
!
! Purpose:
!   - used to obtain the subauroral drift [v-north,v-east,v-vertical]
!     from SAPS electric field model
!
! Calling Sequence:
!   call subaur_drift(umlt, ualat, ukp, vn1, ve1, vv1)
!
! Inputs:
!   umlt  - magnetic local time
!   ualat - auroral latitudes
!   ukp   - Kp index
!
! Outputs:
!   vn1   - SAPS northward ion drift
!   ve1   - SAPS eastward  ion drift
!   vv1   - SAPS vertical  ion drift

    real, intent(in) :: umlt, ualat, ukp
    real, intent(out) :: vn1, ve1, vv1

    integer :: i, j, k
    real :: &
      yn1_new, ye1_new, yv1_new, &
      yn2_new, ye2_new, yv2_new, &
      yn3_new, ye3_new, yv3_new
    real, dimension(nmlt)  :: y1, yn1, ye1, yv1
    real, dimension(nalat) :: y2, yn2, ye2, yv2
    real, dimension(nkp)   :: y3, yn3, ye3, yv3

    do i = 1, nmlt
      do j = 1, nalat
        do k = 1, nkp
          yn3(k) = vn(i,j,k)
          ye3(k) = ve(i,j,k)
          yv3(k) = vv(i,j,k)
        enddo

        call spline(kp,yn3,nkp,y3)
        call splint(kp,yn3,y3,nkp,ukp,yn3_new)

        call spline(kp,ye3,nkp,y3)
        call splint(kp,ye3,y3,nkp,ukp,ye3_new)

        call spline(kp,yv3,nkp,y3)
        call splint(kp,yv3,y3,nkp,ukp,yv3_new)

        yn2(j) = yn3_new
        ye2(j) = ye3_new
        yv2(j) = yv3_new
      enddo

      call spline(alat,yn2,nalat,y2)
      call splint(alat,yn2,y2,nalat,ualat,yn2_new)

      call spline(alat,ye2,nalat,y2) 
      call splint(alat,ye2,y2,nalat,ualat,ye2_new)

      call spline(alat,yv2,nalat,y2)
      call splint(alat,yv2,y2,nalat,ualat,yv2_new)

      yn1(i) = yn2_new
      ye1(i) = ye2_new
      yv1(i) = yv2_new
    enddo

    call spline(mlt,yn1,nmlt,y1)
    call splint(mlt,yn1,y1,nmlt,umlt,yn1_new)

    call spline(mlt,ye1,nmlt,y1)
    call splint(mlt,ye1,y1,nmlt,umlt,ye1_new)

    call spline(mlt,yv1,nmlt,y1)
    call splint(mlt,yv1,y1,nmlt,umlt,yv1_new)

    vn1 = yn1_new
    ve1 = ye1_new
    vv1 = yv1_new

  end subroutine subaur_drift
!----------------------------------------------------------------------------
  subroutine spline(x, y, n, y2)
!   subroutine spline
!   - used for cubic spline interpolation
!   - adopted from Numerical Recipes

    integer, intent(in) :: n
    real, intent(in), dimension(n) :: x, y
    real, intent(out), dimension(n) :: y2

    integer :: i
    real :: p, qn, sig, un
    real, dimension(n) :: u

    y2(1) = 0.
    u(1)  = 0.

    do i=2,n-1
      sig   = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
      p     = sig * y2(i-1) + 2.
      y2(i) = (sig-1.) / p
      u(i)  = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo

    qn = 0.
    un = 0.

    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do i = n-1,1,-1
      y2(i) = y2(i)*y2(i+1)+u(i)
    enddo

  end subroutine spline
!-----------------------------------------------------------------------
  subroutine splint(xa, ya, y2a, n, x, y)
!   subroutine splint
!   - used to call subroutine spline for interpolations
!   - adopted from Numerical Recipes

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
  subroutine handle_error(funcname, ncerr)

    use netcdf, only: nf90_strerror

    character(len=*), intent(in) :: funcname
    integer, intent(in) :: ncerr

    write(6, "('NetCDF error encountered: ', a, ', when calling ', a)") &
      trim(nf90_strerror(ncerr)), funcname

  end subroutine handle_error
!-----------------------------------------------------------------------
end module subaur_module
