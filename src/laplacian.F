!
      subroutine laplacian(lap_arg,lap_out,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Compute Laplacian using Spectral Approximation.
!
      use params_module,only: nlon,nlonp2,nlonp4,nlat
      use addfld_module,only: addfld
      use sh_coef_module,only: lpmn,zmn
      use filter_module,only: ntrigs,trigs,ifax
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
      integer :: i,lat,m,n,nx,nw
      integer,parameter :: nmax=35
      real,dimension(0:nmax-1,0:nmax) :: amn ! a(m,n) spectral coefficient
      real,dimension(nmax-1,nmax)     :: bmn ! b(m,n) spectral coefficient
      real :: fx(nlonp4,nlat),               ! Used by FFT subroutine
     |        wfft((nlonp4+1)*nlat)
!      real :: temp(lon0:lon1,lat0:lat1)      ! used to debug
!
! FFT
! Load fx from f for the fft:
!      temp(:,:) = 0.
      fx(:,:) = 0.
!
!      temp(lon0+2:lon1-2,lat0:lat1) = fx(1:nlon,1:nlat)
!      temp(lon0+2,lat0) = 1.
!      temp(lon0+2,lat1) = 1.
!      temp(lon1-2,lat0) = 1.
!      temp(lon1-2,lat1) = 1.
!      call addfld('TEST1',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      fx(1:nlon,1:nlat) = lap_arg(3:nlonp2,1:nlat)
      !do lat=1,nlat
      !  do i=1,nlon
      !    fx(i,lat) = lap_arg(i,lat)
      !  enddo ! i
      !enddo ! lat
!
!      temp(lon0+2:lon1-2,lat0:lat1) = fx(1:nlon,1:nlat)
!      call addfld('TEST2',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
      nx = (nlonp4)*nlat
      nw = (nlonp4+1)*nlat ! fft999 claims this should be (nlon+1)*nlat
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlat,
     |  -1)
!
!      temp(lon0+2:lon1-2,lat0:lat1) = fx(1:nlon,1:nlat)
!      call addfld('TEST3',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Fit 
      amn(:,:) = 0.
      bmn(:,:) = 0.
      ! n=0,nmax-1
      do n=0,nmax-1
        ! m=0
        do lat=lat0,lat1
          amn(0,n) = amn(0,n)+zmn(lat,0,n)*fx(1,lat)
        enddo ! lat=lat0,lat1
        ! m=1,n
        do m=1,n
          do lat=lat0,lat1
            amn(m,n) = amn(m,n)+zmn(lat,m,n)*fx(2*m+1,lat)
            bmn(m,n) = bmn(m,n)+zmn(lat,m,n)*fx(2*m+2,lat)
          enddo ! lat=lat0,lat1
        enddo ! m=1,n
      enddo ! n=0,nmax-1
      ! n=nmax,m=0
      do lat=lat0,lat1
        amn(0,nmax) = amn(0,nmax)+zmn(lat,0,nmax)*fx(1,lat)
      enddo ! lat=lat0,lat1
      ! n=nmax,m=2,nmax-1,2
      do m=2,nmax-1,2
        do lat=lat0,lat1
          amn(m,nmax) = amn(m,nmax)+zmn(lat,m,nmax)*fx(2*m+1,lat)
          bmn(m,nmax) = bmn(m,nmax)+zmn(lat,m,nmax)*fx(2*m+2,lat)
        enddo ! lat=lat0,lat1
      enddo ! m=2,nmax-1,2
!
!      temp(:,:) = 0.
!      temp(lon0+3:lon0+nmax+1,lat0+1:lat0+nmax) = bmn(1:nmax-1,1:nmax)
!      call addfld('BMN',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!      temp(:,:) = 0.
!      temp(lon0+2:lon0+nmax+1,lat0:lat0+nmax) = amn(0:nmax-1,0:nmax)
!      call addfld('AMN',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Synthesis
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
      fx(:,:) = 0.
      ! n=0,nmax-1
      do n=0,8 !12 !nmax-1
        ! m=0
        do lat=lat0,lat1
          fx(1,lat) = fx(1,lat)+amn(0,n)*lpmn(lat,0,n) ! A(0)
        enddo ! lat=lat0,lat1
        ! m=1,n
        do m=1,n
          do lat=lat0,lat1
            fx(2*m+1,lat) = fx(2*m+1,lat)+amn(m,n)*lpmn(lat,m,n) ! A(1),A(2),...
            fx(2*m+2,lat) = fx(2*m+2,lat)+bmn(m,n)*lpmn(lat,m,n) ! B(1),B(2),...
          enddo ! lat=lat0,lat1
        enddo ! m=1,n
      enddo ! n=0,nmax-1
      ! n=nmax,m=0
      !do lat=lat0,lat1
      !  fx(1,lat) = fx(1,lat)+amn(0,nmax)*lpmn(lat,0,nmax) ! A(0)
      !enddo ! lat=lat0,lat1
      ! n=nmax,m=2,nmax-1,2
      !do m=2,nmax-1,2
      !  do lat=lat0,lat1
      !    fx(2*m+1,lat) = fx(2*m+1,lat)+amn(m,nmax)*lpmn(lat,m,nmax) ! A(2),A(4),...
      !    fx(2*m+2,lat) = fx(2*m+2,lat)+bmn(m,nmax)*lpmn(lat,m,nmax) ! B(2),B(4),...
      !  enddo ! lat=lat0,lat1
      !enddo ! m=2,nmax-1,2
!
!      temp(lon0+2:lon1-2,lat0:lat1) = fx(1:nlon,1:nlat)
!      call addfld('TEST4',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Inverse transform fourier back to gridpoint:
      call fftrans(fx,nx,wfft,nw,trigs,ntrigs,ifax,1,nlonp4,nlon,nlat,
     |  1)
!
!      temp(lon0+2:lon1-2,lat0:lat1) = fx(1:nlon,1:nlat)
!      call addfld('TEST5',' ',' ',
!     |  temp(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Return lap_out from fx:
      lap_out(:,:) = 0.
      lap_out(3:nlonp2,1:nlat) = fx(1:nlon,1:nlat)
!
      end subroutine laplacian
!-----------------------------------------------------------------------
