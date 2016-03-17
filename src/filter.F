!
      module filter_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlon,nlonp1,nlonp2,nlonp4,nlat
      implicit none
!
! Coefficients and factors for fft.
! Setfft is called once per run from init (see init_module.F)
!   Setfft calls set99, which returns trigs and ifax. 
!   Trigs and ifax are used by fftrans (see util.F) calls by filter.f 
!     and filter2.f.
!
      integer,parameter :: ntrigs = 3*nlon/2+1
      real :: trigs(ntrigs) ! e.g., if nlon==72, ntrigs==109
      integer :: ifax(13)
!
      contains
!-----------------------------------------------------------------------
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
      integer :: n1,n2,k,kk,i,ii,nlevs,nx,nw
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
      nx = nlonp4*(2*(lev1-lev0+1))
      nw = (nlonp4+1)*2*(lev1-lev0+1)
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  -1)
!
! Remove wave numbers greater than kutj
      do k = 1,nlevs
        do i=n1,n2
          fx(i,k) = 0.0
        enddo
      enddo
!
! Inverse transform fourier back to gridpoint:
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  1)
!
! Redefine f from fx:
      do k=lev0,lev1
        do i=1,nlon
          f(i+2,k) = fx(i,k)
        enddo
      enddo
      end subroutine filter
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine filter2(f,lev0,lev1,kutj,lat)
      use cons_module,only: dlamda
      use params_module,only: nlon,dlon
!
! Args:
      integer,intent(in) :: kutj,lev0,lev1,lat
      real,intent(inout) :: f(nlonp4,lev0:lev1)
!
! Local:
      integer :: nn(nlat),nlevs,k,i,j,nx,nw,n1,kut2(nlat)
      real :: fx(nlonp4,2*(lev1-lev0+1)),
     |        wfft((nlonp4+1)*2*(lev1-lev0+1))
      real :: smoothfunc ! smoothing function
      real :: coslon
!
! 1/8/08 btf: Define nn for 5.0 or 2.5 degree resolution:
!
#if (NLAT==36 && NLON==72)
      nn=(/43,18, 9, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 9,18,43/) ! 5.0 deg
#elif (NLAT==72 && NLON==144)
      nn=(/200,90,40,40,22,22,14,14,10,10, 8, 8, 6, 6, 4, 4, 2, 2,
     |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     |      2, 2, 4, 4, 6, 6, 8, 8,10,10,14,14,22,22,40,40,90,200/) ! 2.5 deg
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif

#if (NLAT==36 && NLON==72)
      kut2=(/ 0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0, 
     |        0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0, 
     |        0, 0,  0,  0,  0,  0/) ! 5.0 deg
#elif (NLAT==72 && NLON==144)
      kut2=(/0, 0, 1, 2, 4, 4, 6, 6, 8, 8,10,10,12,12,15,15,18,18,
     |      20,20,20,20,18,18,15,12, 8, 8, 4, 4, 4, 4, 2, 2, 1, 1,
     |       1, 1, 2, 2, 4, 4, 4, 4, 8, 8,12,15,18,18,20,20,20,20,
     |      18,18,15,15,12,12,10,10, 8, 8, 6, 6, 4, 4, 2, 1, 0, 0/) ! 2.5 deg
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
! Load field data into fx (1-nlon <- 3-nlonp2)
      nlevs = lev1-lev0+1
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
! subroutine fftrans(a,na,work,nw,trigs,ntrigs,ifax,inc,jump,n,lot,isign)
!
      nx = nlonp4*(2*(lev1-lev0+1))
      nw = (nlonp4+1)*2*(lev1-lev0+1)
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  -1)
!
! Change filters so that it does not over filtering at high latitudes, it will be the
! same as filter for low wavenumer, but wrapping up smoothly for large wavenumers, not a
! sharp transition, so there is still filtering effect in the lower latitudes
!                                                       Wenbin Wang  06/11/13 
      n1=2*(kut2(lat)-1)+3
!
! Multiply by smoothing function:
! Test coslon to avoid underflow in smoothfunc at the poles
!
      do k=lev0,lev1
        do i=n1,nlon
          coslon = cos(((i+1-n1)/2.)*dlamda/2.)
          if (coslon >= 0.1) then
            smoothfunc = coslon**(2*nn(lat))
            fx(i+1,k) = fx(i+1,k)*smoothfunc
          else
            fx(i+1,k) = 0.
          endif
        enddo ! i=1,nlon
      enddo ! k=lev0,lev1
!
! Inverse transform back to gridpoint:
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlevs,
     |  +1)
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
