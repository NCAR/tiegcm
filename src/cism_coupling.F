#if defined(INTERCOMM) || defined(CISMAH)
!-----------------------------------------------------------------------
      module cism_coupling_module
!-----------------------------------------------------------------------
!DESCRIPTION:
!
! Module for coupling TIEGCM with LFM and MIX (CMIT).  Coupling can use
! either InterComm or disk I/O file exchanges (i.e. "adhoc" mode).
! For implementation details, see the source code for the corresponding
! coupling infrastructure;
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Preprocessor flag  !!  Source Code        !!  Notes               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    CISMAH           !!  cism_adhoc.F       !! slow adhoc/file exchg!!
!!    INTERCOMM        !!  cism_intercomm.F   !! High-performance     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      use addfld_module,only: addfld
      use params_module,only: 
     |     nlat,   ! number of geographic latitudes (at 5 deg, nlat==36)
     |     nlatp1, ! nlat+1
     |     nlon,   ! number of geographic longitudes (at 5 deg, nlon=72)
     |     nlonp1, ! nlon+1
     |     nlonp4  ! nlon+4
      integer,parameter :: mytid=0,ntask=1
      real, dimension(nlat,nlonp1) ::
     |     ting_pot_interp,  ! potential in geo coord from MIX
     |     ting_eng_interp,  ! energy in geo coord from MIX
     |     ting_flux_interp  ! flux in geo coord from MIX
      real, dimension(nlat,nlonp1) ::
     |     ting_ped,     ! Temp array --> gzigm1
     |     ting_hall,    ! Temp array --> gzigm2
     |     ting_gnsrhs  ! Temp array --> gnsrhs
      real, dimension(nlonp4,nlat) ::
     |     gpot, ! potential in geographic coordinates, periodic bound.
     |     geng, ! energy in geographic coordinates, periodic boundary
     |     gflx  ! flux in geographic coordinates, periodic boundary
      real, dimension(nlonp1,0:nlat+1) ::
     |     gpotm 

#if defined(INTERCOMM)
      integer :: xjd
#endif
      contains 
!-----------------------------------------------------------------------
      subroutine initialize
#if defined(INTERCOMM)
      call ci_init
#elif defined(CISMAH)
      call ca_init
#endif
      end subroutine initialize
!-----------------------------------------------------------------------
      subroutine import
      integer ii,jj

      write(6,*)'inside receive'
#if defined(INTERCOMM)
      call ci_receive
#elif defined(CISMAH)
      call ca_receive
#endif
!
! Process the imported data
!

!
! Add period points
!
      do jj=1,nlat
       do ii=1,nlonp1-1
        gpot(ii+2,jj)=ting_pot_interp(jj,ii)  !/1.2
        geng(ii+2,jj)=ting_eng_interp(jj,ii)
        gflx(ii+2,jj)=ting_flux_interp(jj,ii)
       enddo
       do ii=1,nlonp1
        gpotm(ii,jj)=ting_pot_interp(jj,ii)
       enddo
       do ii=1,2
        gpot(ii,jj)=gpot(nlonp1+ii-1,jj)
        geng(ii,jj)=geng(nlonp1+ii-1,jj)
        gflx(ii,jj)=gflx(nlonp1+ii-1,jj)
        gpot(nlonp1+ii+1,jj)=gpot(ii,jj)
        geng(nlonp1+ii+1,jj)=geng(ii,jj)
        gflx(nlonp1+ii+1,jj)=gflx(ii,jj)
       enddo
      enddo
!
! Add pole values
!
      polev1=0.
      polev2=0.
      do i=1,nlonp1
       polev1=polev1+gpotm(i,1)
       polev2=polev2+gpotm(i,2)
      enddo
!
! Extend in longitude
!
      do ii=1,nlonp1
       gpotm(ii,0)=(9.*polev1-polev2)/(8.*float(nlonp1))
      enddo

      polev1=0.
      polev2=0.
      do i=1,nlonp1
       polev1=polev1+gpotm(i,nlat)
       polev2=polev2+gpotm(i,nlat-1)
      enddo
      do ii=1,nlonp1
       gpotm(ii,nlat+1)=(9.*polev1-polev2)/(8.*float(nlonp1))
      enddo

      end subroutine import
!-----------------------------------------------------------------------
      subroutine export(modeltime)
! ... Shared Module Variables ..........................................
      use dynamo_module,only: gzigm2,gzigm1,gnsrhs 
! ... Local variables ..................................................
      integer :: jj,ii
! ... Parameter variables ..............................................
      integer, intent(in) :: modeltime(4)      
! ... Begin ............................................................
!
! Prepare data for export:
!     
      do jj=1,nlat
         ! Add periodic points
         gzigm1(nlonp1,jj)=gzigm1(1,jj)
         gzigm2(nlonp1,jj)=gzigm2(1,jj)
         gnsrhs(nlonp1,jj)=gnsrhs(1,jj)         
         
         do ii=1,nlonp1
            if(gzigm1(ii,jj)<0.2)gzigm1(ii,jj)=0.2 !wwb
            if(gzigm2(ii,jj)<0.2)gzigm2(ii,jj)=0.2 !wwb      
            
!           FIXME:  Copying data into temporary arrays:
!                   ting_ped,ting_hall,ting_gnsrhs ... and then copying to 
!                   the arrays gzigm1,gzigm2,gnsrhs in cism_coupling.F
!                   This code should be refactored to prevent ting_* tmp arrays.
            ting_ped(jj,ii) = gzigm1(ii,jj)
            ting_hall(jj,ii) = gzigm2(ii,jj)
            ting_gnsrhs(jj,ii) = gnsrhs(ii,jj)
         enddo
      enddo

      write(6,*) "Sending at: ", 
     1     modeltime(1),modeltime(2), modeltime(3),modeltime(4)
      
!
! Export the data...
!      
#if defined(INTERCOMM)
      call ci_send(modeltime)
#elif defined(CISMAH)
      call ca_send(modeltime)
#endif
!
! Save exchange variables to secondary history files.
!
      call cism_save
      end subroutine export
!-----------------------------------------------------------------------
      subroutine finalize
#if defined(INTERCOMM)
      call ci_close
#elif defined(CISMAH)
      call ca_close
#endif
      end subroutine finalize
!-----------------------------------------------------------------------
#if 0
! FIXME:  This code is un-tested with CMIT 2.5.  It's designed to 
!         scatter CMIT exchange data to all the nodes in an MPI run of
!         the TIEGCM.
      subroutine cism_scatter(f,i0,i1,j0,j1)
! 
! Redistribute data from task 1 to other tasks
!
! Uses:
!
      use params_module,only: nlonp4,nlat
      use mpi_module,only: mxlon,mxlat,mytidi,mytidj,ntaski,ntaskj,
     |  itask_table,mytid,tasks,mytid
#ifdef MPI
#include <mpif.h>
      integer ::
     |  irstat(MPI_STATUS_SIZE)      ! mpi receive status
#endif
#ifndef MPI
      integer :: irstat,mpi_comm_world,mpi_real8
#endif
!
! Args:
!
      integer,intent(in) :: i0,i1,j0,j1
      real,intent(inout) :: f(nlonp4,nlat)
!
! Local:
!
      integer :: i,j,nlonsend,nlatsend,len,idest,isrc,ier,istart,
     |  isend,irecv,itask,jtask,lonsend0,lonsend1,latsend0,latsend1,mtag
      real,allocatable ::
     |  sndbuf(:,:), ! send buffer (mxlon,mxlat)
     |  rcvbuf(:,:)  ! recv buffer (mxlon,mxlat)
!      implicit none
!
! Allocate send and receive buffers:
!
      allocate(sndbuf(mxlon,mxlat),stat=ier)
      if (ier /= 0) 
     |  write(6,"('>>> mp_scatterlons_f3d: error allocating sndbuf.')")
      allocate(rcvbuf(mxlon,mxlat),stat=ier)
      if (ier /= 0) 
     |  write(6,"('>>> mp_scatterlons_f3d: error allocating rcvbuf.')")
      len = mxlon*mxlat
!
! Send data from mytidi==0 and mytidj==1 to other processes
!
      if(mytid==0)then
       do jtask=0,ntaskj-1
        if (jtask==0)then
         istart=1
        else
         istart=0
        endif
        do itask=istart,ntaski-1
         idest = itask_table(itask,jtask)
         lonsend0 = tasks(idest)%lon0
         lonsend1 = tasks(idest)%lon1
         nlonsend = lonsend1-lonsend0+1
         latsend0 = tasks(idest)%lat0
         latsend1 = tasks(idest)%lat1
         nlatsend = latsend1-latsend0+1
         sndbuf(1:nlonsend,1:nlatsend) = 
     |       f(lonsend0:lonsend1,latsend0:latsend1) 
         mtag = idest+mytid
         call mpi_isend(sndbuf,len,MPI_REAL8,idest,mtag,MPI_COMM_WORLD,
     |    isend,ier)
c        if (ier /= 0) 
c     |    call handle_mpi_err(ier,'mp_scatterlons_f3d send to idest')
         call mpi_wait(isend,irstat,ier)
c        if (ier /= 0) 
c     |    call handle_mpi_err(ier,'mp_scatterlons_f3d wait for send')
        enddo ! itask=istart,ntaski-1
       enddo ! jtask=0,ntaskj-1
      else
!
! Receive data from task (0,0):
!
       isrc = itask_table(0,0)
       mtag = isrc+mytid
       call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,mtag,MPI_COMM_WORLD,
     |    irecv,ier)
c      if (ier /= 0) 
c     |  call handle_mpi_err(ier,'mp_scatterlons_f3d recv fm isrc')
       call mpi_wait(irecv,irstat,ier)
c      if (ier /= 0) 
c      |  call handle_mpi_err(ier,'mp_scatterlons_f3d wait for recv')
       do j=j0,j1 
         f(i0:i1,j) = rcvbuf(1:i1-i0+1,j-j0+1)
       enddo ! j=j0,j1
      endif
!
! Free local buffer space:
!
      deallocate(sndbuf)
      deallocate(rcvbuf)
      
      end subroutine cism_scatter
! end #if 0 ...
! FIXME: Implement & test cism_scatter for TIEGCM MPI coupling.
#endif
!-----------------------------------------------------------------------
      subroutine cism_save
! This subroutine saves physical parameters that are either important to 
! the M-I coupling physics or crucial for code debugging.
!
! Add these variables to secondary history file:
!    1. high latitude potential form M-I coupler (2D)
!    2. high latitude precipitation characteristic energy from M-I
!       coupler (2D)
!    3. high latitude precipitation particle number flux from M-I 
!       coupler (2D)
!    4. global height-integrated Pedersen conductance from dynamo.F 
!       (2D,73,0:37)
!    5. global height-integrated Hall conductanbce from dynamo.F
!       (2D,73,0:37)
!-----------------------------------------------------------------------
      use params_module,only: nlonp1,nlon
      use dynamo_module,only: gzigm2,gzigm1,gnsrhs
      use fields_module,only: ped,hall,  levd0,levd1
!-----------------------------------------------------------------------
 
      call addfld('gpot', 'Potential from M-I Coupler',' ', 
     |     gpot, 'lon',1,nlonp4,  'lat',1,nlat,  0)
      
      call addfld('gpotm', 'Potential from M-I Coupler (mag)',' ', 
     |     gpotm, 'lon',1,nlonp1,  'lat',1,nlat,  0)
      
      call addfld("geng", 'Energy from M-I Coupler', ' ', 
     |     geng, 'lon',1,nlonp4,  'lat',1,nlat,  0)
      
      call addfld("gflx", 'Number Flux from M-I Coupler ', ' ', 
     |     gflx, 'lon',1,nlonp4,  'lat',1,nlat,  0)
      
      call addfld("gzigm1", 'Pedersen Conductance (geographic)', ' ', 
     |     gzigm1, 'lon',1,nlonp1, 'lat',1,nlat,  0)
      
      call addfld("gzigm2", 'Hall Conductance (geographic) ', ' ', 
     |     gzigm2, 'lon',1,nlonp1, 'lat',1,nlat,  0)
      
      call addfld("gnsrhs", 'Height-integrated neutral wind ', ' ', 
     |     gnsrhs, 'lon',1,nlonp1,  'lat',1,nlat,  0)
      
!
! FIXME:  Something is wrong with these ped/hall addfld calls... 
!         Maybe because lat/lon dims are wrong (off by +/- 4)?
!         See r783 of this file for the original implementation.
!
!      call addfld("ped", 'Altitude profile of Pedersen Cond ', ' ', 
!     |     ped, 'lev',levd0,levd1,  'lon',1,nlonp1,  nlat)      
!      
!      call addfld("hall", 'Altitue profile of Hall Cond. ', ' ', 
!     |     hall, 'lev',levd0,levd1,  'lon',1,nlonp1,  nlat)
      end subroutine cism_save
!-----------------------------------------------------------------------
! FIXME:  Should move cism_ucurrent from dynamo.F to here!
!-----------------------------------------------------------------------
      end module cism_coupling_module
!-----------------------------------------------------------------------
      subroutine cism_pot2mag
!
! This subroutine converts high latitude potential in geographic coordinate 
! obtained form the CMIT M-I couplier to geomagnetci coordinate to be used
! in dynamo calculations
!
! Uses
!
      use params_module,only: nlonp1,nmlonp1,nmlon,nmlat
      use cism_coupling_module, only:gpotm
      use dynamo_module,only: phihm,geo2mag
      use magfield_module,only: ig,jg,wt
!
! Local:
!
      integer :: jj

      do jj=1,nmlat
       call geo2mag
     1   (phihm(1,jj),gpotm,ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,jj)
c       do i=1,nmlon
c        write(6,*)phihm(i,jj),i,jj
c       enddo
      enddo

      end subroutine cism_pot2mag
!-----------------------------------------------------------------------
#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.  
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine cism_coupling_null
      end subroutine cism_coupling_null
!-----------------------------------------------------------------------
#endif
