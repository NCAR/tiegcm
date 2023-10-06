module interp_module

  implicit none

  contains
!-----------------------------------------------------------------------
  pure function interp3d(z,x,y,zp,xp,yp,fp,zlog) result(f)
! 3d interpolation, linear extrapolation
! 3d dimensions in the order of z,x,y in accordance with model fields
! zp,xp,yp: input locations; fp: input values
! z,x,y: output locations; f: output values
! zp,xp,yp are required to be monotonically increasing

    real,dimension(:),intent(in) :: z,x,y,zp,xp,yp
    real,dimension(size(zp),size(xp),size(yp)),intent(in) :: fp
    logical,intent(in),optional :: zlog
    real,dimension(size(z),size(x),size(y)) :: f

    integer :: k,nx,ny,i,j
    real,dimension(size(zp),size(x),size(y)) :: aux

    forall (k=1:size(zp)) aux(k,:,:) = interp2d(x,y,xp,yp,fp(k,:,:))

    nx = size(x)
    ny = size(y)

! logarithmic interpolation in z direction if zlog is set to true
    if (present(zlog)) then
      if (zlog) then
        forall (i=1:nx,j=1:ny) f(:,i,j) = exp(interp1d(z,zp,log(aux(:,i,j))))
        return
      endif
    endif

! otherwise (zlog not present or set to false), linear interpolation
    forall (i=1:nx,j=1:ny) f(:,i,j) = interp1d(z,zp,aux(:,i,j))

  end function interp3d
!-----------------------------------------------------------------------
  pure function interp2d(x,y,xp,yp,fp) result(f)
! 2d interpolation, linear extrapolation
! xp,yp: input locations; fp: input values
! x,y: output locations; f: output values
! xp,yp are required to be monotonically increasing

    real,dimension(:),intent(in) :: x,y,xp,yp
    real,dimension(size(xp),size(yp)),intent(in) :: fp
    real,dimension(size(x),size(y)) :: f

    integer :: i,j
    real,dimension(size(xp),size(y)) :: aux

    forall (i=1:size(xp)) aux(i,:) = interp1d(y,yp,fp(i,:))
    forall (j=1:size(y)) f(:,j) = interp1d(x,xp,aux(:,j))

  end function interp2d
!-----------------------------------------------------------------------
  pure function interp1d(x,xp,fp) result(f)
! 1d interpolation, linear extrapolation
! xp: input locations; fp: input values
! x: output locations; f: output values
! xp is required to be monotonically increasing

    real,dimension(:),intent(in) :: x,xp
    real,dimension(size(xp)),intent(in) :: fp
    real,dimension(size(x)) :: f

    integer :: nx,nxp,i
    integer,dimension(size(x)) :: i0,i1
    real,dimension(size(x)) :: x0,x1

    nx = size(x)
    nxp = size(xp)

    forall (i=1:nx) i0(i) = find(xp,x(i))
    i1 = i0+1

    where (i0 <= 0)
      i0 = 1
      i1 = 2
    endwhere
    where (i1 >= nxp+1)
      i0 = nxp-1
      i1 = nxp
    endwhere

    x0 = xp(i0)
    x1 = xp(i1)

    f = ((x1-x)*fp(i0) + (x-x0)*fp(i1)) / (x1-x0)

  end function interp1d
!-----------------------------------------------------------------------
  pure function find(x,x0) result(i)
! 1d binary search
! x: value array; x0: value to be found
! i: the index satisfying x(i)<=x0<x(i+1)
! x is required to be monotonically increasing

    real,dimension(:),intent(in) :: x
    real,intent(in) :: x0
    integer :: i

    integer :: nx,i0,i1

    nx = size(x)

    if (x0 < x(1)) then
      i = 0
    elseif (x0 >= x(nx)) then
      i = nx
    else
      i0 = 1
      i1 = nx
      do while (i0+1 < i1)
        i = (i0+i1)/2
        if (x(i) <= x0) then
          i0 = i
        else
          i1 = i
        endif
      enddo
      i = (i0+i1)/2
    endif

  end function find
!-----------------------------------------------------------------------
end module interp_module
