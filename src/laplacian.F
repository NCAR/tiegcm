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
      use params_module,only: nlon,nlonp2,nlonp4,nlat,glon,glat
      use addfld_module,only: addfld
      use mpi_module,only: mp_gather_f2d,mp_scatter_f2d,tasks,ntask,
     |  mytid
!
! Select correct set of coefficients for He calculation, depending
! on grid resolution:
!
      use he_coefs_module,only: nlat_he,nlon_he,glat_he,glon_he,
     |  trigs,ntrigs,ifax,pmn,zmn
      use interp_module,only: interp2d
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
      integer :: j,lat,m,n,nx,nw,
     |  truncdeg ! Degree of truncation to maintain numerical stability
      integer,parameter :: nmax=35
      real,dimension(0:nmax-1,0:nmax) :: amn ! a(m,n) spectral coefficient
      real,dimension(nmax-1,nmax)     :: bmn ! b(m,n) spectral coefficient
      real :: fx(nlonp4,nlat),               ! Used by FFT subroutine
     |        wfft((nlon_he+5)*nlat_he)
      real :: fx_tmp(nlonp4,nlat)
      real,dimension(nlon_he+4,nlat_he) :: fx_he
      real,dimension(nlon+1) :: glonp1
      real,dimension(nlon_he+1) :: glonp1_he

      if (nlat==36 .and. nlon==72) then
        truncdeg = 4
      else
        truncdeg = 8
      endif

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
      call mp_gather_f2d(lap_arg,fx_tmp,tasks(:)%lon0,tasks(:)%lon1,
     |  tasks(:)%lat0,tasks(:)%lat1,ntask,nlonp4,nlat,1)
      fx(1:nlon,1:nlat) = fx_tmp(3:nlonp2,1:nlat)
!
! Root task does the fft:
!
      if (mytid==0) then
!
        if (nlat == nlat_he) then
          fx_he = fx
        else
          glonp1(1:nlon) = glon
          glonp1(nlon+1) = glon(1)+360
          fx(nlon+1,:) = fx(1,:)
          fx_he(1:nlon_he,:) =
     |      interp2d(glon_he,glat_he,glonp1,glat,fx(1:nlon+1,:))
          fx_he(nlon_he+1:nlon_he+4,:) = 0.
        endif
!
! Forward transform gridpoint to fourier:
! (fftrans is in util.F)
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
        nx = (nlon_he+4)*nlat_he
        nw = (nlon_he+5)*nlat_he ! fft999 claims this should be (nlon+1)*nlat
        call fftrans(fx_he,nx,wfft,nw,trigs,ntrigs,ifax,1,nlon_he+4,
     |    nlon_he,nlat_he,-1)
!
! Fit (now on global grid)
        amn(:,:) = 0.
        bmn(:,:) = 0.
        ! n=0,nmax-1
        do n=0,truncdeg ! for no truncation, sum to nmax-1
          ! m=0
          do lat=1,nlat_he
            amn(0,n) = amn(0,n)+zmn(lat,0,n)*fx_he(1,lat)
          enddo 
          ! m=1,n
          do m=1,n
            do lat=1,nlat_he
              amn(m,n) = amn(m,n)+zmn(lat,m,n)*fx_he(2*m+1,lat)
              bmn(m,n) = bmn(m,n)+zmn(lat,m,n)*fx_he(2*m+2,lat)
            enddo
          enddo ! m=1,n
        enddo ! n=0,nmax-1
        ! for (no) truncation (un)comment the following loops
        ! n=nmax;m=0
!        do lat=1,nlat_he
!          amn(0,nmax) = amn(0,nmax)+zmn(lat,0,nmax)*fx_he(1,lat)
!        enddo
        ! n=nmax;m=2,nmax-1,2
!        do m=2,nmax-1,2
!          do lat=1,nlat_he
!            amn(m,nmax) = amn(m,nmax)+zmn(lat,m,nmax)*fx_he(2*m+1,lat)
!            bmn(m,nmax) = bmn(m,nmax)+zmn(lat,m,nmax)*fx_he(2*m+2,lat)
!          enddo
!        enddo ! m=2,nmax-1,2
!
! Synthesis
!   ORDERING OF COEFFICIENTS:
!     A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!     WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
        fx_he(:,:) = 0.
        ! n=0,nmax-1
        do n=0,truncdeg ! for no truncation, sum to nmax-1
          ! m=0
          do lat=1,nlat_he
            fx_he(1,lat) = fx_he(1,lat)-n*(n+1)*amn(0,n)*pmn(lat,0,n)   ! A(0)
          enddo
          ! m=1,n
          do m=1,n
            do lat=1,nlat_he
              fx_he(2*m+1,lat) = fx_he(2*m+1,lat)-n*(n+1)*amn(m,n)*
     |                        pmn(lat,m,n)                              ! A(1),A(2),...
              fx_he(2*m+2,lat) = fx_he(2*m+2,lat)-n*(n+1)*bmn(m,n)*
     |                        pmn(lat,m,n)                              ! B(1),B(2),...
            enddo
          enddo ! m=1,n
        enddo ! n=0,nmax-1
        ! for (no) truncation (un)comment the following loops
        ! n=nmax;m=0
!        do lat=lat0,lat1
!          fx_he(1,lat) = fx_he(1,lat)-nmax*(nmax+1)*amn(0,nmax)*
!     |                pmn(lat,0,nmax)                                   ! A(0)
!        enddo ! lat=lat0,lat1
        ! n=nmax;m=2,nmax-1,2
!        do m=2,nmax-1,2
!          do lat=lat0,lat1
!            fx_he(2*m+1,lat) = fx_he(2*m+1,lat)-nmax*(nmax+1)*amn(m,nmax)*
!     |                      pmn(lat,m,nmax)                             ! A(2),A(4),...
!            fx_he(2*m+2,lat) = fx_he(2*m+2,lat)-nmax*(nmax+1)*bmn(m,nmax)*
!     |                      pmn(lat,m,nmax)                             ! B(2),B(4),...
!          enddo ! lat=lat0,lat1
!        enddo ! m=2,nmax-1,2
!
! Inverse transform fourier back to gridpoint:
        call fftrans(fx_he,nx,wfft,nw,trigs,ntrigs,ifax,1,nlon_he+4,
     |    nlon_he,nlat_he,1)

        if (nlat == nlat_he) then
          fx = fx_he
        else
          glonp1_he(1:nlon_he) = glon_he
          glonp1_he(nlon_he+1) = glon_he(1)+360
          fx_he(nlon_he+1,:) = fx_he(1,:)
          fx(1:nlon,:) =
     |      interp2d(glon,glat,glonp1_he,glat_he,fx_he(1:nlon_he+1,:))
        endif

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
!
! Scatter results to subdomains:
!
      call mp_scatter_f2d(fx_tmp,lap_out,tasks(:)%lon0,tasks(:)%lon1,
     |  tasks(:)%lat0,tasks(:)%lat1,ntask,nlonp4,nlat,1)

!     do j=lat0,lat1
!       write(6,"('laplacian: j=',i4,' lap_out(:,j)=',/,(6e12.4))")
!    |    j,lap_out(:,j)
!     enddo

!     call addfld('LAP_OUT',' ',' ',
!    |  lap_out,'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      end subroutine laplacian
!-----------------------------------------------------------------------
