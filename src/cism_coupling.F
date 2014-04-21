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
      use aurora_module,only:	! (nlonp1,0:nlatp1)
     |     gekev, ! electron energy in keV in geographic coordinates (also at poles for mag2geo)
     |     gmwpm2 ! electron energy flux mW/m2 in geographic coordinates (also at poles for mag2geo)

      real*8, dimension(nlat,nlonp1) ::
     |     ting_pot_interp,  ! potential in geo coord from MIX
     |     ting_eng_interp,  ! energy in geo coord from MIX
     |     ting_flux_interp  ! flux in geo coord from MIX
      real*8, dimension(nlat,nlonp1) ::
     |     ting_ped,     ! Temp array --> gzigm1
     |     ting_hall,    ! Temp array --> gzigm2
     |     ting_gnsrhs  ! Temp array --> gnsrhs
      real*8, dimension(nlonp4,nlat) ::
     |     gpot ! potential in geographic coordinates, periodic bound.
!     |     geng, ! characteristic energy in geographic coordinates, periodic boundary
!     |     gflx  ! number flux in geographic coordinates, periodic boundary

!     Note: gpotm latitude dimension is defined to match mag2geo
!     specifications.
      real*8, dimension(nlonp4,0:nlatp1) ::
     |     gpotm ! Potential in geographic coord from M-I coupler

#if defined(INTERCOMM)
      integer :: xjd
#endif
      contains 
!-----------------------------------------------------------------------
      subroutine initialize
      use mpi_module,only: mytid,ntask
#ifndef MPI
      mytid = 0
      ntask = 1
#endif

#if defined(INTERCOMM)
      call ci_init
#elif defined(CISMAH)
      call ca_init
#endif
      end subroutine initialize
!-----------------------------------------------------------------------
      subroutine import
!!      subroutine import (model_hpi)
! 1.  Receive arrays
!       - ting_pot_interp
!       - ting_eng_interp
!       - ting_flx_interp
!     from MIX on process zero.
! 2.  Broadcast to all MPI tasks
! 3.  Set to gpot, gpotm, geng, gflx (gpot,m with periodic boundary; 
!     geng,gflx,gpotm at poles also)
! 4.  Clean up pole values for Dynamo solver.
!
      use input_module,only: ctpoten
      use mpi_module,only: mytid
#ifdef MPI
      use mpi_module,only: mytid, handle_mpi_err
! Args:
#include <mpif.h>
      integer iErr
#endif
!!      real,intent(out) :: model_hpi(2)
      integer ii,jj,js,jn
      real cpmaxsh,cpminsh,cpmaxnh,cpminnh
! 11/08 EMERY added eflux=flux*fac_p2e*alfa   (fac_p2e=2*1.602e-9, so 2*alfa in this factor)
      real, parameter:: fac_p2e = 3.204e-9 ! convert from particle to energy flux

      write(6,"(A,I1,A)") "TIEGCM (",mytid,"): inside receive"

#if defined(INTERCOMM)
      call ci_receive
#elif defined(CISMAH)
      call ca_receive
#endif

#ifdef MPI
!     FIXME: could optimize this by packing all three arrays
!     (ting_pot_interp, ting_eng_interp, ting_flux_interp) into one 3d
!     array and calling MPI_Bcast with count=3*nlat*nlonp1.
      call MPI_Bcast(ting_pot_interp, nlat*nlonp1, MPI_REAL8, 0, 
     $     MPI_COMM_WORLD, iErr)
      if (iErr .ne. 0) then
         call handle_mpi_err(iErr, 'cism_coupling.F: call mpi_bcast')
      endif

      call MPI_Bcast(ting_eng_interp, nlat*nlonp1, MPI_REAL8, 0, 
     $     MPI_COMM_WORLD, iErr)
      if (iErr .ne. 0) then
         call handle_mpi_err(iErr, 'cism_coupling.F: call mpi_bcast')
      endif

      call MPI_Bcast(ting_flux_interp, nlat*nlonp1, MPI_REAL8, 0, 
     $     MPI_COMM_WORLD, iErr)
      if (iErr .ne. 0) then
         call handle_mpi_err(iErr, 'cism_coupling.F: call mpi_bcast')
      endif
#endif

!
! Process the imported data.
!

      do jj=1,nlat
         do ii=1,nlon
            ! Note ii+2 skips over period points, where geng, and gflx don't have them
! NOTE:  ting_eng_interp is divided by 2 from MIX before entering TIEGCM so geng is alfa 
!        ting_flux_interp is the number flux from MIX so gflx is electron number flux
            gpot(ii+2,jj)=ting_pot_interp(jj,ii) !/1.2
!           geng(ii+2,jj)=ting_eng_interp(jj,ii)
!           gflx(ii+2,jj)=ting_flux_interp(jj,ii)
         enddo
         do ii=1,nlonp1
            gpotm(ii+2,jj)=ting_pot_interp(jj,ii)
! NOTE:  ting_eng_interp is divided by 2 from MIX before entering TIEGCM so multiply 
!         gekev by 2 to get electron energy in keV 
!        ting_flux_interp is the number flux from MIX so convert gmwpm2 to energy flux 
!         in mW/m2 as gmwpm2=ting_flux_interp*fac_p2e*gekev/2 (where gekev/2=ting_eng_interp)
            gekev(ii,jj)=ting_eng_interp(jj,ii) * 2.
            gmwpm2(ii,jj)=ting_flux_interp(jj,ii)*fac_p2e*
     |        ting_eng_interp(jj,ii)
         enddo
      enddo	

!     
! Dynamo solver requires pole values for gpotm.
! Linearly interpolate to add pole values
!
! Southern Hemisphere:  jj=latitudie=0
!
      polev1=0.
      polev2=0.
      do i=1,nlonp1
         polev1=polev1+gpotm(i+2,1)
         polev2=polev2+gpotm(i+2,2)
      enddo

      do ii=1,nlonp1
         gpotm(ii+2,0)=(9.*polev1-polev2)/(8.*float(nlonp1))
      enddo

!
! Northern Hemisphere: jj=nLat+1
!
      polev1=0.
      polev2=0.
      do i=1,nlonp1
         polev1=polev1+gpotm(i+2,nlat)
         polev2=polev2+gpotm(i+2,nlat-1)
      enddo
      do ii=1,nlonp1
         gpotm(ii+2,nlatp1)=(9.*polev1-polev2)/(8.*float(nlonp1))
      enddo

!
! Set periodic points
!     
      do jj=1,nlat
         do ii=1,2
            gpot(ii,jj)=gpot(nlon+ii,jj)
!            geng(ii,jj)=geng(nlon+ii,jj)
!            gflx(ii,jj)=gflx(nlon+ii,jj)
            gpotm(ii,jj)=gpotm(nlon+ii,jj)
            
            gpot(nlonp1+ii+1,jj)=gpot(ii+2,jj)
!            geng(nlonp1+ii+1,jj)=geng(ii+2,jj)
!            gflx(nlonp1+ii+1,jj)=gflx(ii+2,jj)
            gpotm(nlonp1+ii+1,jj)=gpotm(ii+2,jj)
         enddo
      enddo

! 08/13 bae:  Added call to calculate HP from geflux in SH and NH
! 03/14 bae:  Deleted since model_hpi now found in aurora.F from gflx
!     call cism_input_analysis (model_hpi)

!!!
!!!  08/11:  FIXME: While merging TIEGCM 1.94.1 into LTR-2.1.4-beta, I found
!!!          the following code segment.  Does CTPOTEN really need to be
!!!          calculated here?  What does CTPOTEN do in the input namelist?
!!!          Isn't CTPOTEN read from the GPI file?  Is CMIT ready for the 
!!!          dynamic crit mods?  Commenting this code out for now.  Need
!!!          to talk with Ben on this.  Compare to TIEGCM revisions r571
!!!          and r575 (TIEGCM repository).
!!!
!!!!
!!!! 01/11:  Find ctpoten (kV) or min/max average from both hemispheres (gpot in V)
!!!!         ctpoten is NOT (yet) used to find theta0 in aurora_cons, where theta0
!!!!         is used in colath for crit(1,2).  Set theta0=10 so crit1,2=15,30 deg (old).
!!!!
!!!      cpmaxsh = -1000000.
!!!      cpmaxnh = -1000000.
!!!      cpminsh =  1000000.
!!!      cpminnh =  1000000.
!!!
!!!      do js=1,nlat/2
!!!       jn=nlat/2+js
!!!       do ii=1,nlonp1+1
!!!	cpmaxsh = max(cpmaxsh,gpot(ii,js))
!!!	cpminsh = min(cpminsh,gpot(ii,js))
!!!	cpmaxnh = max(cpmaxnh,gpot(ii,jn))
!!!	cpminnh = min(cpminnh,gpot(ii,jn))
!!!       enddo
!!!      enddo
!!!      ctpoten = 0.5*(cpmaxsh-cpminsh+cpmaxnh-cpminnh)*0.001
!!!!     write (6,"(1x,'cism CP (SH,NH,av) =',3f8.2)") 
!!!!   |  (cpmaxsh-cpminsh)*0.001,(cpmaxnh-cpminnh)*0.001,ctpoten	

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
         do ii=1,2
            gzigm1(ii,jj)=gzigm1(nlon+ii,jj)
            gzigm2(ii,jj)=gzigm2(nlon+ii,jj)
            gnsrhs(ii,jj)=gnsrhs(nlon+ii,jj)
            gzigm1(nlonp1+ii+1,jj)=gzigm1(ii,jj)
            gzigm2(nlonp1+ii+1,jj)=gzigm2(ii,jj)
            gnsrhs(nlonp1+ii+1,jj)=gnsrhs(ii,jj)
         enddo
         
         do ii=1,nlonp1
            if(gzigm1(ii+2,jj)<0.2)gzigm1(ii+2,jj)=0.2 !wwb
            if(gzigm2(ii+2,jj)<0.2)gzigm2(ii+2,jj)=0.2 !wwb      
            
!           FIXME:  Copying data into temporary arrays:
!                   ting_ped,ting_hall,ting_gnsrhs ... and then copying to 
!                   the arrays gzigm1,gzigm2,gnsrhs in cism_coupling.F
!                   This code should be refactored to prevent ting_* tmp arrays.
            ting_ped(jj,ii) = gzigm1(ii+2,jj)
            ting_hall(jj,ii) = gzigm2(ii+2,jj)
            ting_gnsrhs(jj,ii) = gnsrhs(ii+2,jj)
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
 
      call addfld('gpot','Potential from M-I Coupler (geographic)',' ', 
     |     gpot, 'lon',1,nlonp4,  'lat',1,nlat,  0)
      
!      call addfld("geng",'Energy from M-I Coupler', ' ',
!     |     geng, 'lon',1,nlonp4,  'lat',1,nlat,  0)
!      
!      call addfld("gflx",'Number Flux from M-I Coupler ', ' ', 
!     |     gflx, 'lon',1,nlonp4,  'lat',1,nlat,  0)

      call addfld("geng",'Energy from M-I Coupler', ' ',
     |     ting_eng_interp, 'lon',1,nlonp1,  'lat',1,nlat,  0)
      
      call addfld("gflx",'Number Flux from M-I Coupler ', ' ', 
     |     ting_flux_interp, 'lon',1,nlonp1,  'lat',1,nlat,  0)

      call addfld("gekev",'Energy*2 from M-I Coupler', ' ',
     |     gekev, 'lon',1,nlonp1,  'lat',1,nlat,  0)
      
      call addfld("gmwpm2",'Energy Flux from M-I Coupler ', ' ', 
     |     gmwpm2, 'lon',1,nlonp1,  'lat',1,nlat,  0)

      call addfld('gpotm','Potential from M-I Coupler (geographic)',' ', 
     |     gpotm(:, 1:nlat), 'lon',1,nlonp4,  'lat',1,nlat,  0)

      call addfld("gzigm1",'Pedersen Conductance (geographic)', ' ', 
     |     gzigm1(:, 1:nlat), 'lon',1,nlonp4, 'lat',1,nlat,  0)
      
      call addfld("gzigm2",'Hall Conductance (geographic) ', ' ', 
     |     gzigm2(:, 1:nlat), 'lon',1,nlonp4, 'lat',1,nlat,  0)

      call addfld("gnsrhs",'Height-integrated neutral wind ', ' ', 
     |     gnsrhs(:, 1:nlat), 'lon',1,nlonp4, 'lat',1,nlat,  0)
      
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
         call geo2mag(phihm(1,jj),gpotm(3:nlonp1+2,:),
     |        ig,jg,wt,nlonp1,nmlonp1,nmlon,nmlat,jj)
         ! Periodic point
         phihm(nmlonp1,jj) = phihm(1,jj)
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
