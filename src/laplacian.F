#include <defs.h>
      subroutine laplacian(lap_arg,lap_out,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Compute scalar surface Laplacian:  Nabla^2 F = 1/sin(theta)*d/d(theta)*
!  (sin(theta)*dF/d(theta)) + 1/sin^2(theta)*d^2*F/d(phi)^2 using a 
!  truncated spectral approximation, where theta = colatitude, phi = 
!  longitude, and F includes the constant factor 1/re^2 that is normally 
!  present on the rhs.
!
      use params_module,only: nlon,nlonp2,nlonp4,nlat
      use addfld_module,only: addfld
      use filter_module,only: ntrigs,trigs,ifax
#ifdef MPI
      use mpi_module,only: mp_gather_f2d,mp_scatter_f2d,tasks,ntask,
     |  mytid
#endif
!
! Select correct set of coefficients for He calculation, depending
! on grid resolution:
!
#if (NLAT==36 && NLON==72)
      use he_coefs_sres,only: pmn,zmn  ! he_coefs_sres.F
#elif (NLAT==72 && NLON==144)
      use he_coef0_dres,only: pmn      ! he_coef0_dres.F
      use he_coef1_dres,only: zmn      ! he_coef1_dres.F
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,dimension(lon0:lon1,lat0:lat1),intent(in)::
     |  lap_arg    ! Laplacian argument
      real,dimension(lon0:lon1,lat0:lat1),intent(out)::
     |  lap_out    ! Laplacian output
!
! Local:
#ifndef MPI
      integer,save :: mytid=0
#endif
#if (NLAT==36 && NLON==72)
      integer,parameter :: truncdeg=4 ! Degree of truncation to maintain numerical stability
#elif (NLAT==72 && NLON==144)
      integer,parameter :: truncdeg=8 ! Degree of truncation to maintain numerical stability
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
      integer :: i,j,lat,m,n,nx,nw,lonbeg,lonend
      integer,parameter :: nmax=35
      real,dimension(0:nmax-1,0:nmax) :: amn ! a(m,n) spectral coefficient
      real,dimension(nmax-1,nmax)     :: bmn ! b(m,n) spectral coefficient
      real :: fx(nlonp4,nlat),               ! Used by FFT subroutine
     |        wfft((nlonp4+1)*nlat)
      real :: fx_tmp(nlonp4,nlat)
      real :: lap_tmp(lon0:lon1,lat0:lat1)

!     call addfld('LAP_ARG',' ',' ',lap_arg,'lon',lon0,lon1,
!    |  'lat',lat0,lat1,0)
!
! FFT
! Load fx from f for the fft:
      fx(:,:) = 0.
!
! Gather 2d subdomain arrays to the root task. 
! (fx_tmp is returned global array on root task)
!
#ifdef MPI
      call mp_gather_f2d(lap_arg,fx_tmp,tasks(:)%lon0,tasks(:)%lon1,
     |  tasks(:)%lat0,tasks(:)%lat1,ntask,nlonp4,nlat,1)
      fx(1:nlon,1:nlat) = fx_tmp(3:nlonp2,1:nlat)
#else
      fx(1:nlon,1:nlat) = lap_arg(3:nlonp2,1:nlat)
#endif
!
! Root task does the fft:
!
      if (mytid==0) then
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
        nx = (nlonp4)*nlat
        nw = (nlonp4+1)*nlat ! fft999 claims this should be (nlon+1)*nlat
        call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlat,
     |    -1)
!
! Fit (now on global grid)
        amn(:,:) = 0.
        bmn(:,:) = 0.
        ! n=0,nmax-1
        do n=0,truncdeg ! for no truncation, sum to nmax-1
          ! m=0
          do lat=1,nlat
            amn(0,n) = amn(0,n)+zmn(lat,0,n)*fx(1,lat)
          enddo 
          ! m=1,n
          do m=1,n
            do lat=1,nlat
              amn(m,n) = amn(m,n)+zmn(lat,m,n)*fx(2*m+1,lat)
              bmn(m,n) = bmn(m,n)+zmn(lat,m,n)*fx(2*m+2,lat)
            enddo
          enddo ! m=1,n
        enddo ! n=0,nmax-1
        ! for (no) truncation (un)comment the following loops
        ! n=nmax;m=0
!        do lat=1,nlat
!          amn(0,nmax) = amn(0,nmax)+zmn(lat,0,nmax)*fx(1,lat)
!        enddo
        ! n=nmax;m=2,nmax-1,2
!        do m=2,nmax-1,2
!          do lat=1,nlat
!            amn(m,nmax) = amn(m,nmax)+zmn(lat,m,nmax)*fx(2*m+1,lat)
!            bmn(m,nmax) = bmn(m,nmax)+zmn(lat,m,nmax)*fx(2*m+2,lat)
!          enddo
!        enddo ! m=2,nmax-1,2
!
! Synthesis
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
        fx(:,:) = 0.
        ! n=0,nmax-1
        do n=0,truncdeg ! for no truncation, sum to nmax-1
          ! m=0
          do lat=1,nlat
            fx(1,lat) = fx(1,lat)-n*(n+1)*amn(0,n)*pmn(lat,0,n)         ! A(0)
          enddo
          ! m=1,n
          do m=1,n
            do lat=1,nlat
              fx(2*m+1,lat) = fx(2*m+1,lat)-n*(n+1)*amn(m,n)*
     |                        pmn(lat,m,n)                              ! A(1),A(2),...
              fx(2*m+2,lat) = fx(2*m+2,lat)-n*(n+1)*bmn(m,n)*
     |                        pmn(lat,m,n)                              ! B(1),B(2),...
            enddo
          enddo ! m=1,n
        enddo ! n=0,nmax-1
        ! for (no) truncation (un)comment the following loops
        ! n=nmax;m=0
!        do lat=lat0,lat1
!          fx(1,lat) = fx(1,lat)-nmax*(nmax+1)*amn(0,nmax)*
!     |                pmn(lat,0,nmax)                                   ! A(0)
!        enddo ! lat=lat0,lat1
        ! n=nmax;m=2,nmax-1,2
!        do m=2,nmax-1,2
!          do lat=lat0,lat1
!            fx(2*m+1,lat) = fx(2*m+1,lat)-nmax*(nmax+1)*amn(m,nmax)*
!     |                      pmn(lat,m,nmax)                             ! A(2),A(4),...
!            fx(2*m+2,lat) = fx(2*m+2,lat)-nmax*(nmax+1)*bmn(m,nmax)*
!     |                      pmn(lat,m,nmax)                             ! B(2),B(4),...
!          enddo ! lat=lat0,lat1
!        enddo ! m=2,nmax-1,2
!
! Inverse transform fourier back to gridpoint:
        call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlat,
     |    1)

        fx_tmp = 0.
        do j=1,nlat
          fx_tmp(3:nlonp2,j) = fx(1:nlon,j)
        enddo

!       call addfld('LAP_FX',' ',' ',fx_tmp,'lon',1,nlonp4,
!    |    'lat',1,nlat,0)
!
! End root task only:
      endif ! mytid==0
!
! Return lap_out from fx:
!
      lap_out(:,:) = 0. ! init
#ifdef MPI
!
! Scatter results to subdomains:
!
      call mp_scatter_f2d(fx_tmp,lap_out,tasks(:)%lon0,tasks(:)%lon1,
     |  tasks(:)%lat0,tasks(:)%lat1,ntask,nlonp4,nlat,1)

#else
      lap_out(3:nlonp2,1:nlat) = fx(1:nlon,1:nlat)
#endif
!     do j=lat0,lat1
!       write(6,"('laplacian: j=',i4,' lap_out(:,j)=',/,(6e12.4))")
!    |    j,lap_out(:,j)
!     enddo

!     call addfld('LAP_OUT',' ',' ',
!    |  lap_out,'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      end subroutine laplacian
!-----------------------------------------------------------------------
