!
      module filter_module
      use params_module,only: nlon,nlonp1,nlonp2,nlonp4,nlat
      implicit none
!
! Coefficients and factors for fft.
! Setfft is called once per run from init (see init_module.F)
!   Setfft calls set99, which returns trigs and ifix. 
!   Trigs and ifax are used by fftrans (see util.F) calls by filter.f 
!     and filter2.f.
!
      real :: trigs(3*nlon/2+1)
      integer :: ifax(13)
!
      contains
!-----------------------------------------------------------------------
!
      subroutine filter(f,lev0,lev1,kutj,lat)
!
! Remove longitudinal waves of prognostic variables with global fft.
! Remove wave numbers greater than kutj (see kut(nlat) in cons_module.F)
! This is called after mp_gatherlons, and only by tasks with mytidi==0.
! On entry, task must have global longitude data defined (mp_gatherlons).
!
! Args:
      integer,intent(in) :: kutj,lev0,lev1,lat
      real,intent(inout) :: f(nlonp4,lev0:lev1)
!
! Local:
      integer :: n1,n2,k,kk,i,ii,nlevs
      real :: fx(nlonp4,2*(lev1-lev0+1)),
     |        wfft((nlonp4+1)*2*(lev1-lev0+1))
!
      nlevs = lev1-lev0+1
      n1 = 2*kutj+3 ! nyquist freq (?)
      n2 = nlon+2
      if (n1 > n2) then
!       write(6,"('filter: lat=',i2,' kutj=',i2,' n1,2=',
!    |    2i3,' n1 > n2')") lat,kutj,n1,n2
        return
      endif
!
! Load fx from f for the fft:
      fx(:,:) = 0.
      do k=lev0,lev1
        do i=1,nlon
          fx(i,k) = f(i+2,k)
        enddo
      enddo
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
      call fftrans(fx,wfft,trigs,ifax,1,nlonp4,nlon,nlevs,-1)
!
! Remove wave numbers greater than kutj
      do k = 1,nlevs
        do i=n1,n2
          fx(i,k) = 0.0
        enddo
      enddo
!
! Inverse transform fourier back to gridpoint:
      call fftrans(fx,wfft,trigs,ifax,1,nlonp4,nlon,nlevs,1)
!
! Redefine f from fx:
      do k=lev0,lev1
        do i=1,nlon
          f(i+2,k) = fx(i,k)
        enddo
      enddo
      end subroutine filter
!-----------------------------------------------------------------------
      subroutine filter2(f,lev0,lev1,kutj,lat)
      use cons_module,only: dlamda
      use params_module,only: nlon,dlon
!
! Args:
      integer,intent(in) :: kutj,lev0,lev1,lat
      real,intent(inout) :: f(nlonp4,lev0:lev1)
!
! Local:
      integer,parameter :: nlond4 = nlon/4
      integer :: nn(nlond4),nlevs,k,i
      integer :: nn_5(72/4) ! for 5 deg in lon
      real :: fx(nlonp4,2*(lev1-lev0+1)),
     |        wfft((nlonp4+1)*2*(lev1-lev0+1))
!
! For 5 degree resolution: 
      nn_5=(/45,20,11, 7, 5,4,3,2,2,1,1,1,1,1,1,1,1,1/) ! 5 deg grid
      if (dlon==5.0) then
        do i=1,nlond4
          nn(i)=nn_5(i)
        enddo
      elseif (dlon==2.5) then ! 2.5 deg res
        do i=1,nlond4
          nn(i*2-1) = nn_5(i)
          nn(i*2) = nn_5(i)
        enddo
!       write(6,"('filter2: nlon=',i3,' nlon/4=',i3,' nn=',/,
!    |    (12i4))") nlon,nlon/4,nn
      else
        write(6,"('filter2: unsupported dlon=',f8.2)") dlon
        stop 'dlon'
      endif

!
      nlevs = lev1-lev0+1
      if (kutj > nlat/2) return
      do k=lev0,lev1
        do i=1,nlon
          fx(i,k) = f(i+2,k)
        enddo
        do i=nlonp1,nlonp2
          fx(i,k)=0.
        enddo
      enddo
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
      call fftrans(fx,wfft,trigs,ifax,1,nlonp4,nlon,nlevs,-1)
!
! Multiply by smoothing function:
      do k=lev0,lev1
        do i=1,nlon
          fx(i+2,k) = fx(i+2,k)*cos(((i+1)/2)*dlamda/2.)**(2*nn(kutj))
        enddo ! i=1,nlon
      enddo ! k=lev0,lev1
!
! Inverse transform back to gridpoint:
      call fftrans(fx,wfft,trigs,ifax,1,nlonp4,nlon,nlevs,1)
!
! Save smoothed field:
      do k=lev0,lev1
        do i=1,nlon
          f(i+2,k) = fx(i,k)
        enddo ! i=1,nlon
      enddo ! k=lev0,lev1
      end subroutine filter2
!-----------------------------------------------------------------------
      end module filter_module
