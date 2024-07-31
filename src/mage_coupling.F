!-----------------------------------------------------------------------
      module mage_coupling_module
!-----------------------------------------------------------------------
! DESCRIPTION:

! Module for MPI coupling TIEGCM with GAMERA and REMIX (MAGE).

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Preprocessor flag  !!  Source Code        !!  Notes               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GAMERA             !!  mage_coupling.F    !!  MPI                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      use params_module,only: nlonp4,nlat,nlatp1,
     |  nmlatp1,nmlonp1,nmlat,nmlon,rp
      use input_module,only: oneway
      !use pdynamo_module,only: phihm
#ifdef GAMERA
      use mage_oneway
#endif
      implicit none

      real(rp),dimension(nlonp4,nlat) ::
     |  geng,        ! energy in geographic coordinates, periodic boundary
     |  gflx,        ! flux in geographic coordinates, periodic boundary
     |  gpot         ! potential in geographic coordinates, periodic boundary
!    |  gBbeAvgEng,  ! BBE Precip Average Eng-geo, periodic boundary
!    |  gBbeNumFlux, ! BBE Precip Num flux-geo, periodic boundary
!    |  gCuspAvgEng, ! Cusp Precip Average Eng-geo, periodic boundary
!    |  gCuspNumFlux ! Cusp Precip Num flux-geo, periodic boundary

! Potential in geographic coordinates from M-I coupler
! Note: gpotm latitude dimension is defined to match mag2geo specifications
      real(rp),dimension(nlonp4,0:nlatp1) :: gpotm
! Values on APEX grid
      real(rp),dimension(nmlonp1,nmlat) :: apotm,aflxm,
     |  grad_x,grad_y,grad_mag

      real(rp) :: aurllbN,aurllbS ! all are co-lats in radians

! Change the following parameters if adding more variables
      integer,parameter ::
     |  nmixinapex = 2, ! APEX Potential, Flux (for auroral bc)
     |  nmixingeo  = 2, ! GEO Flux, Energy
     |  nmixoutapex= 2, ! APEX Pedersen, Hall
     |  nmixoutgeo = 0, ! No GEO exports
     |  nhoutvar   = 9, ! TN, UN, VN, OMEGA, O2, O1, NO, Z, HE
     |  nhinvar    = 0  !

      real(rp),dimension(nmlonp1,2) :: aurorabc1d

      contains 
!-----------------------------------------------------------------------
      subroutine initialize_coupling

! This subroutine initialize MPI interface for data exchange
! Send REMIX the grid information

      use params_module,only: nlatp2,nlonp1,glat,glon,
     |  nmlatp1,nmlonp1,gmlat,gmlon
      use cons_module,only: dtr
      use mpi_module,only: CplCommSize,CplRank,IAm,myAppId,
     |  voltId,gamId,rcmId,hidraNId,hidraSId,hidraId,
     |  hidraNCplRank,hidraSCplRank,hidraCplRank,
     |  mytid

      integer :: ierr,i
      real(rp),dimension(nlatp2) :: glatp2
      real(rp),dimension(nlonp1) :: glonp1

      !gpot = 0
      geng = 0.5_rp
      gflx = 0.0_rp
      gpotm = 0.0_rp
      apotm = 0.0_rp
      aflxm  = 0.0_rp
      !aurllbN = 10*dtr
      !aurllbS = 10*dtr

! Hard code the poles
      glatp2(2:nlatp2-1) = glat
      glatp2(1) = -90.0_rp
      glatp2(nlatp2) = 90.0_rp

! Periodic boundary
      glonp1(1:nlonp1-1) = glon
      glonp1(nlonp1) = glon(1)+360.0_rp

#if defined (GAMERA) || defined (HIDRA)
      do i=1,CplCommSize
            ! Skip Self
            if (i == CplRank+1) continue
            ! Assign rank if match
            select case (IAm(i))
              case (voltId)
                call init_mpi_remix(glatp2,glonp1,gmlat,gmlon)
              case (gamId)
                !write(*,*) "T not coupling to Gam yet"
              case (rcmId)
                !write(*,*) "T not coupling to RCM yet"
              case (hidraNId)
                !write(*,*) "T COUPLING TO HIDRAN"
                call init_mpi_hidra(myAppId+hidraNId,hidraNCplRank)
              case (hidraSId)
                !write(*,*) "T COUPLING TO HIDRAS"
                call init_mpi_hidra(myAppId+hidraSId,hidraSCplRank)
              case (hidraId)
                !write(*,*) "T COUPLING TO HIDRA"
                call init_mpi_hidra(myAppId+hidraId,hidraCplRank)
              case (myAppId)
                !write(*,*) "T is T"
              case default
                if (IAm(i) .eq. 0) cycle
                if (mytid == 0)
     |  write(*,*) "T does not know about this Coupling ID: ",
     |              IAm(i)
              end select
      enddo
      write(*,*) "T Done Coupling Init"
#endif

      end subroutine initialize_coupling
!-----------------------------------------------------------------------
      subroutine import_mage

! 1. Receive arrays (pot, eng, flx) from REMIX on process zero
! 2. Broadcast to all MPI tasks
! 3. Set periodic boundaries for gpot, geng, gflx, gpotm
! 4. Clean up pole values for the dynamo solver

      use params_module,only: nlatp2,nlon,nlonp1,nlonp2
      use mpi_module,only: mytid,TIEGCM_WORLD,handle_mpi_err
      use mpi

      integer :: ierr,i,j
      real(rp),dimension(nlatp2,nlonp1,nmixingeo) :: gvar2d
      real(rp),dimension(nmlat,nmlonp1,nmixinapex) :: avar2d

!     write(6,"('TIEGCM (',i3,'): inside receive')") mytid

      avar2d = 0.0_rp
      gvar2d = 0.0_rp

#ifdef GAMERA
      if (.not. oneway) then
      if (mytid == 0) then

        call import_remix(avar2d,gvar2d)
! Unit conversion
        avar2d(:,:,1) = avar2d(:,:,1)*1e3_rp ! potential: kV -> V
        gvar2d(:,:,1) = gvar2d(:,:,1)*0.5_rp ! mean energy -> characteristic energy
      endif

      call mpi_bcast(avar2d, nmlat*nmlonp1*nmixinapex,
     |  MPI_REAL8, 0, TIEGCM_WORLD, ierr)
      if (ierr /= 0) call handle_mpi_err(ierr, 'import: TIEGCM bcast')

      call mpi_bcast(gvar2d, nlatp2*nlonp1*nmixingeo,
     |  MPI_REAL8, 0, TIEGCM_WORLD, ierr)
      if (ierr /= 0) call handle_mpi_err(ierr, 'import: TIEGCM bcast')

!      call mpi_bcast(aurorabc1d, nmlonp1*2,
!     |  MPI_REAL8, 0, TIEGCM_WORLD, ierr)
!      if (ierr /= 0) call handle_mpi_err(ierr, 'import: TIEGCM bcast')

! Process the imported data
      do j = 1,nlat
        do i = 1,nlon
          !gpot(i+2,j) = var2d(j+1,i,1)
          geng(i+2,j) = gvar2d(j+1,i,1)
          gflx(i+2,j) = gvar2d(j+1,i,2)
          !gpotm(i+2,j) = var2d(j+1,i,1)
        enddo
        !gpotm(nlonp1+2,j) = var2d(j+1,nlonp1,1)
      enddo
      do j=1,nmlat
        do i=1,nmlonp1
          ! nmlonp1,nmlat
          apotm(i,j) = avar2d(j,i,1)
          aflxm(i,j) = avar2d(j,i,2)
        enddo
      enddo
      else
        call ca_receive(gpot,geng,gflx,gpotm)
      endif
#endif
      end subroutine import_mage
!-----------------------------------------------------------------------

      subroutine process_import

#ifdef GAMERA
      use params_module,only: nlatp2,nlon,nlonp1,nlonp2,gmlat,gmlon
      use input_module,only: ctpoten
      use cons_module,only: pi,dlatm,dlonm
      use mpi_module,only: mytid,handle_mpi_err,TIEGCM_WORLD
      use mpi

      integer :: v,ierr,i,j,js,jn
      real(rp) :: polev1,polev2,
     |   cpmaxsh,cpminsh,cpmaxnh,cpminnh,
     |   temp
      real(rp) :: nfluxllb = 1.0e7_rp
      real(rp) :: potllb = 5000._rp ! 1 kV
      real(rp),dimension(nlon) :: aurllbj
  

! Dynamo solver requires pole values for gpotm
! Linearly interpolate to add pole values

! Southern Hemisphere: j=0
      !polev1=0.
      !polev2=0.
      !do i=1,nlonp1
      !  polev1=polev1+gpotm(i+2,1)
      !  polev2=polev2+gpotm(i+2,2)
      !enddo

      !do i=1,nlonp1
      !  gpotm(i+2,0)=(9.*polev1-polev2)/(8.*float(nlonp1))
      !enddo

! Northern Hemisphere: j=nlatp1
      !polev1=0.
      !polev2=0.
      !do i=1,nlonp1
      !  polev1=polev1+gpotm(i+2,nlat)
      !  polev2=polev2+gpotm(i+2,nlat-1)
      !enddo
      !do i=1,nlonp1
      !  gpotm(i+2,nlatp1)=(9.*polev1-polev2)/(8.*float(nlonp1))
      !enddo

! Set periodic points
      do j = 1,nlat
        do i = 1,2
          !gpot(i,j) = gpot(nlon+i,j)
          geng(i,j) = geng(nlon+i,j)
          gflx(i,j) = gflx(nlon+i,j)
          !gpotm(i,j) = gpotm(nlon+i,j)

          !gpot(nlonp2+i,j) = gpot(i+2,j)
          geng(nlonp2+i,j) = geng(i+2,j)
          gflx(nlonp2+i,j) = gflx(i+2,j)
          !gpotm(nlonp2+i,j) = gpotm(i+2,j)
        enddo
      enddo

!!!  08/11:  FIXME: While merging TIEGCM 1.94.1 into LTR-2.1.4-beta, I found
!!!          the following code segment.  Does CTPOTEN really need to be
!!!          calculated here?  What does CTPOTEN do in the input namelist?
!!!          Isn't CTPOTEN read from the GPI file?  Is MAGE ready for the 
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

!!!      cpmaxsh = maxval(gpot(:,1:nlat/2))
!!!      cpmaxnh = maxval(gpot(:,nlat/2+1:nlat))
!!!      cpminsh = minval(gpot(:,1:nlat/2))
!!!      cpminnh = minval(gpot(:,nlat/2+1:nlat))

!!!      ctpoten = 0.5*(cpmaxsh-cpminsh+cpmaxnh-cpminnh)*0.001

      ! Find the lowest lat of all MLT where mixFlux drops to below certain threshold.
      ! Pass it into colath.F, replace crit(1) and 2.
      ! Make NH/SH separately.
      ! Note in cism_adhoc.F, allocate(mixPot(nmixlat,nmixlon)), allocate(mixlon(nmixlon,nmixlat))

      ! real*8,dimension(nlonp1,2) :: aurorabc1d

      aurllbN = 0.0_rp
      aurllbS = 0.0_rp
      aurorabc1d = 0.0_rp
      do j=1,nmlonp1   ! NORTH
        do i=nmlat/2+1,nmlat
          if(aflxm(j,i)>=nfluxllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
        enddo
        i = min(i,nmlat)
        aurorabc1d(j,1) = max(90.0_rp-gmlat(i),15._rp) ! aurllbj is co-lat.
      enddo
!!!      aurllbN = maxval(aurllbj)

      do j=1,nmlonp1   ! SOUTH
        do i=nmlat/2,1,-1
          if(aflxm(j,i)>=nfluxllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
        enddo
        i = max(i,1)
        aurorabc1d(j,2) = max(90._rp+gmlat(i),15._rp) ! aurllbj is co-lat from south pole! Backwards.
      enddo
!!!      aurllbS = maxval(aurllbj)
        aurllbN = maxval(aurorabc1d(:,1))
        aurllbS = maxval(aurorabc1d(:,2))
        !write(*,*) "TGCMCPL: aurora: ",aurllbN,aurllbS

!        ! Let's try a new method
!        aurorabc1d = 0.0_rp
!        grad_x = 0.0_rp
!        grad_y = 0.0_rp
!        grad_mag = 0.0_rp
!
!        ! dumb gradient
!        do j = 2,nmlon
!          do i=nmlat/2+2,nmlat-1
!            grad_x(j,i) = (apotm(j+1,i) - apotm(j-1,i))/(2.0_rp)!*dlonm)
!            grad_y(j,i) = (apotm(j,i+1) - apotm(j,i-1))/(2.0_rp)!*dlatm)
!          end do
!          do i=nmlat/2-1,2,-1
!            grad_x(j,i) = (apotm(j+1,i) - apotm(j-1,i))/(2.0_rp)!*dlonm)
!            grad_y(j,i) = (apotm(j,i-1) - apotm(j,i+1))/(2.0_rp)!*dlatm)
!          end do
!        end do
!
!        ! Compute the magnitude of the gradient
!        do j = 2, nmlon
!            do i = nmlat/2+2,nmlat-1
!                grad_mag(i,j) = sqrt(grad_x(j,i)**2 + grad_y(j,i)**2)
!            end do
!        end do
!
!        ! Calculate gradient boundary
!        do j=1,nmlonp1   ! NORTH
!          do i=nmlat/2+1,nmlat
!            if(grad_mag(j,i)>=potllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
!          enddo
!          i = min(i,nmlat)
!          aurorabc1d(j,1) = max(90.0_rp-gmlat(i),15._rp) ! aurllbj is co-lat.
!        enddo
!!!!      aurllbN = maxval(aurllbj)
!
!        do j=1,nmlonp1   ! SOUTH
!          do i=nmlat/2,1,-1
!            if(grad_mag(j,i)>=potllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
!          enddo
!          i = max(i,1)
!          aurorabc1d(j,2) = max(90._rp+gmlat(i),15._rp) ! aurllbj is co-lat from south pole! Backwards.
!        enddo
!
!        ! Keep the higher latitude between auroral and pot boundary
!        temp = maxval(aurorabc1d(:,1))
!        aurllbN = min(temp,aurllbN)
!        temp = maxval(aurorabc1d(:,2))
!        aurllbS = min(temp,aurllbS)


#endif

      end subroutine process_import

!-----------------------------------------------------------------------
      subroutine export_mage(modeltime)
        use mpi_module,only: CplCommSize,CplRank,IAm,
     |  voltId,gamId,rcmId,hidraNId,hidraSId,hidraId,
     |  hidraNCplRank,hidraSCplRank,hidraCplRank,
     |  myAppId,mytid

      integer,dimension(4),intent(in) :: modeltime
      logical :: hidra_prep

      real(rp),dimension(:,:,:), allocatable :: mixgeoout,mixapexout
      real(rp),dimension(:,:,:,:), allocatable :: hidraout

      integer :: nreq,i
#if defined (GAMERA) || defined (HIDRA)
      ! Prepare the export data
      hidra_prep = .false.
      nreq = 0
      if (mytid == 0) write(*,*) "T Starting Export Prep"
      do i=1,CplCommSize
            ! Skip Self
            if (i == CplRank+1) continue
            ! Assign rank if match
            select case (IAm(i))
              case (voltId)
                call prep_export_remix(modeltime,mixapexout,mixgeoout)
                nreq = nreq + 1
              case (gamId)
                !write(*,*) "T not coupling to Gam yet"
              case (rcmId)
                !write(*,*) "T not coupling to RCM yet"
              case (hidraNId)
                if (.not.hidra_prep) then
                  !write(*,*) "Prep H Export"
                  call prep_export_hidra(modeltime,hidraout)
                  hidra_prep = .true.
                endif
                nreq = nreq + 1
              case (hidraSId)
                if (.not.hidra_prep) then
                  !write(*,*) "Prep H Export"
                  call prep_export_hidra(modeltime,hidraout)
                  hidra_prep = .true.
                endif
                nreq = nreq + 1
              case (hidraId)
                if (.not.hidra_prep) then
                  call prep_export_hidra(modeltime,hidraout)
                  hidra_prep = .true.
                endif
                nreq = nreq + 1
              case (myAppId)
                !write(*,*) "T is T"
              case default
                if (IAm(i) .eq. 0) cycle
                if (mytid == 0)
     |  write(*,*) "T does not know about this Coupling ID: ",
     |              IAm(i)
              end select
      enddo

      ! Allocate request array here


      ! Send the export data
      do i=1,CplCommSize
        ! Skip Self
        if (i == CplRank+1) continue
        ! Assign rank if match
        select case (IAm(i))
          case (voltId)
            call export_remix(mixapexout,mixgeoout)
          case (gamId)
            !write(*,*) "T not coupling to Gam yet"
          case (rcmId)
            !write(*,*) "T not coupling to RCM yet"
          case (hidraNId)
            !write(*,*) "Export to HidraN"
            call export_hidra(hidraout,myAppId+hidraNId,hidraNCplRank)
          case (hidraSId)
            !write(*,*) "Export to HidraS"
            call export_hidra(hidraout,myAppId+hidraSId,hidraSCplRank)
          case (hidraId)
            call export_hidra(hidraout,myAppId+hidraId,hidraCplRank)
          case (myAppId)
            !write(*,*) "T is T"
          case default
                if (IAm(i) .eq. 0) cycle
                if (mytid == 0)
     |  write(*,*) "T does not know about this Coupling ID: ",
     |              IAm(i)
              end select
      enddo

      ! Deallocate arrays
      do i=1,CplCommSize
            ! Skip Self
            if (i == CplRank+1) continue
            ! Assign rank if match
            select case (IAm(i))
              case (voltId)
                if (allocated(mixapexout)) deallocate(mixapexout)
                if (allocated(mixgeoout)) deallocate(mixgeoout)
              case (gamId)
                !write(*,*) "T not coupling to Gam yet"
              case (rcmId)
                !write(*,*) "T not coupling to RCM yet"
              case (hidraNId)
                if (allocated(hidraout)) deallocate(hidraout)
              case (hidraSId)
                if (allocated(hidraout)) deallocate(hidraout)
              case (hidraId)
                if (allocated(hidraout)) deallocate(hidraout)
              case (myAppId)
                !write(*,*) "T is T"
              case default
                if (IAm(i) .eq. 0) cycle
                if (mytid == 0)
     |  write(*,*) "T does not know about this Coupling ID: ",
     |              IAm(i)
              end select
      enddo

      ! Deallocate request array here
#endif
      end subroutine export_mage
!-----------------------------------------------------------------------
      subroutine prep_export_hidra(modeltime,hidra_out)

! ... Shared Module Variables ..........................................
      use params_module,only: nlon,nlonp4,nlat,nlevp1
      use fields_module,only: f4d,itp,itc
      use fields_module,only: i_z,i_tn,i_un,i_vn,i_w,i_o2,i_o1,i_no,
     |  i_he

      use mpi_module,only: lon0,lon1,lat0,lat1,mp_gather_f3d,mytid,
     |  ntask,handle_mpi_err
      use mpi

! ... Parameter variables ..............................................
      integer,dimension(4),intent(in) :: modeltime

! ... Local variables ..................................................
      integer :: i,j,v,k,ierr, ilist(nhoutvar)
      real(rp),dimension(:,:,:,:),allocatable :: f3d_sub,f3d_glb,
     |  hidra_out

! ... Begin ............................................................

      ! Number of fields we are passing
      ! This list must match the order that is in tgcmInterp enumerator
      ilist = (/i_tn, i_un, i_vn, i_he, i_o2, i_o1, i_no, i_w, i_z /)

! Prepare data for export
      !write(*,*) "PREP HIDRA EXPORT STUFF: ",mytid
      if (.not. allocated(f3d_sub))
     |  allocate(f3d_sub(nlevp1,lon0:lon1,lat0:lat1,nhoutvar))
      if (.not. allocated(f3d_glb))
     |  allocate(f3d_glb(nlevp1,nlonp4,nlat,nhoutvar))

      do v =1, nhoutvar
       f3d_sub(:,:,:,v) = f4d(ilist(v))%data(:,lon0:lon1,lat0:lat1,itp)
      enddo
      call mp_gather_f3d(f3d_sub,f3d_glb,
     |  lon0,lon1,lat0,lat1,ntask,nhoutvar)

      if (mytid == 0) then
        write(*,*) "TGCM: ALLOCATE"
        if (.not. allocated(hidra_out))
     |    allocate(hidra_out(nlevp1,nlon,nlat,nhoutvar))

        hidra_out(:,:,:,:) = f3d_glb(:,3:nlon+2,:,:)
        do i=1,nlon
          do j=1,nlat
            do k=1,nlevp1-1
              if ((f3d_glb(k,i,j,1) .lt.1e-5_rp) .or.
     |    (isnan(f3d_glb(k,i,j,1)))) then
              write(*,*) "Bad TN: ",f3d_glb(k,i,j,1),k,i,j
              endif
            enddo
          enddo
        enddo
        ! Check fidelity
        do v=1,nhoutvar
          write(*,*) 'TGCM: max/min V:',maxval(hidra_out(:,:,:,v)),
     |      minval(hidra_out(:,:,:,v))
        enddo

        if (mytid == 0) then
        !write(6,"('TGCM: Sending to HIDRA at: ',4i4)") modeltime
        write(6,"('GCM2HIDRA: ',4i4)") modeltime
        endif
      endif

! Save exchange variables to secondary history files

!     call addfld('gpot','Potential from M-I Coupler (geographic)',' ',
!    |  gpot(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('geng','Energy from M-I Coupler', ' ',
!    |  geng(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('gflx','Number Flux from M-I Coupler', ' ',
!    |  gflx(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('gpotm','Potential from M-I Coupler (geographic)',' ',
!    |  gpotm(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)

      if (allocated(f3d_sub)) deallocate(f3d_sub)
      if (allocated(f3d_glb)) deallocate(f3d_glb)
      end subroutine prep_export_hidra

!-----------------------------------------------------------------------

      subroutine export_hidra(hidra_out,commId,hidraCplRank)
        use params_module,only: nlon,nlonp4,nlat,nlevp1
        use mpi_module,only: CplComm,handle_mpi_err,mytid
        use mpi

        real(rp),dimension(:,:,:,:) :: hidra_out
        integer :: commId,ierr,hidraCplRank

      ! Export the data
        if (CplComm /= MPI_COMM_NULL) then
        write(*,*) "TIEGCM EXPORT TO HIDRA",mytid,commId,hidraCplRank
            call mpi_send(hidra_out, nlevp1*nlon*nlat*nhoutvar,
     |        MPI_DOUBLE_PRECISION, hidraCplRank,
     |        commId*100, CplComm, ierr)
            if (ierr /= 0)
     |        call handle_mpi_err(ierr, 'export: TIEGCM send to REMIX')
            write(*,*) "TGCM: Done sending to HIDRA"
        endif
      end subroutine export_hidra

!-----------------------------------------------------------------------
      subroutine prep_export_remix(modeltime,avar2d,gvar2d)

! ... Shared Module Variables ..........................................
      use params_module,only: nlatp2,nlon,nlonp1,nlonp2
      use fields_module,only: gzigm1,gzigm2,gnsrhs
      use mpi_module,only: mytid,ntask,tasks
      use mpi_module,only: mlat0,mlat1,mlon0,mlon1,
     |  mp_gather_pdyn
      use mpi_module,only: lon0,lon1,lat0,lat1,mp_gather_f2d
      use addfld_module,only: addfld
      use fields_module,only: azigm1,azigm2
      use mpi

! ... Parameter variables ..............................................
      integer,dimension(4),intent(in) :: modeltime

! ... Local variables ..................................................
      integer :: i,j,v,ierr
      real(rp),dimension(:,:,:),allocatable :: avar2d,gvar2d

      real(rp) :: amsub(mlon0:mlon1,mlat0:mlat1,nmixoutapex)
      real(rp) :: amglb(nmlonp1,nmlat,nmixoutapex)

      real(rp),allocatable,dimension(:,:,:) :: gsub
      real(rp),allocatable,dimension(:,:,:) :: gglb

! ... Begin ............................................................

! Prepare data for export
      amsub(:,:,1) = azigm1(mlon0:mlon1,mlat0:mlat1)
      amsub(:,:,2) = azigm2(mlon0:mlon1,mlat0:mlat1)

      where(amsub < 0.2_rp) amsub = 0.2_rp

      call mp_gather_pdyn(amsub,mlon0,mlon1,mlat0,mlat1,
     |  amglb,nmlonp1,nmlat,nmixoutapex)

      if (mytid == 0) then
        if (.not. allocated(avar2d)) then
           allocate(avar2d(nmlat,nmlonp1,nmixoutapex))
        endif
        do v=1,nmixoutapex
          do j = 1,nmlat
            do i = 1,nmlonp1
              avar2d(j,i,v) = amglb(i,j,v)
            enddo
          enddo
          ! Fill pole with average
          !avar2d(1,:,v) = sum(amglb(1:nmlon,1,v))/nmlon
          !avar2d(nmlatp2,:,v) = sum(amglb(1:nmlon,nmlat,v))/nmlon
        enddo
      endif

        if (nmixoutgeo .ne. 0) then
          ! Prepare geographic output arrays
          ! Gather GEO arrays
          if (.not. allocated(gsub))
     |  allocate(gsub(lon0:lon1,lat0:lat1,nmixoutgeo))
          if (.not. allocated(gglb))
     |  allocate(gglb(nlonp4,nlat,nmixoutgeo))
          do v=1,nmixoutgeo
          select case(v)
          case (1)
          gsub(:,:,v) = gzigm1
          case(2)
          gsub(:,:,v) = gzigm2
          where(gsub < 0.2_rp) gsub = 0.2_rp
          case(3)
          gsub(:,:,v) = gnsrhs
          end select
          enddo

          call mp_gather_f2d(gsub(lon0:lon1,lat0:lat1,nmixoutgeo),gglb,
     |  tasks(:)%lon0,tasks(:)%lon1,tasks(:)%lat0,tasks(:)%lat1,
     |  ntask,nlonp4,nlat,nmixoutgeo)

      if (mytid == 0) then
        if (.not. allocated(gvar2d)) then
           allocate(gvar2d(nlatp2,nlonp1,nmixoutgeo))
        endif
          do v=1,nmixoutgeo
            write(*,*) "Preparing GEO stuff here"
            do j = 2,nlatp1
              do i = 1,nlonp1
                gvar2d(j,i,v) = gglb(i+2,j-1,v)
              enddo
            enddo
            gvar2d(1,:,1) = sum(gglb(3:nlonp2,1,v))/nlon
            gvar2d(nlatp2,:,1) = sum(gglb(3:nlonp2,nlat,v))/nlon
          enddo

        endif
      endif

        if (mytid == 0)
     | write(6,"('TGCM: Sending at: ',4i4)") modeltime

! Save exchange variables to secondary history files

!     call addfld('gpot','Potential from M-I Coupler (geographic)',' ',
!    |  gpot(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('geng','Energy from M-I Coupler', ' ',
!    |  geng(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('gflx','Number Flux from M-I Coupler', ' ',
!    |  gflx(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)
!     call addfld('gpotm','Potential from M-I Coupler (geographic)',' ',
!    |  gpotm(lon0:lon1,lat0:lat1),'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine prep_export_remix

!-----------------------------------------------------------------------

      subroutine export_remix(avar2d,gvar2d)
        use params_module,only: nlatp2,nlonp1
        use mpi_module, only: CplComm,mixCplRank,handle_mpi_err,
     |   myAppId,voltId
        use mpi

        real(rp),dimension(:,:,:) :: avar2d,gvar2d
        integer :: ierr

! Export the data
        if (CplComm /= MPI_COMM_NULL) then

          if ( nmixoutgeo .ne. 0) then
          call mpi_send(gvar2d, nlatp2*nlonp1*nmixoutgeo,
     |       MPI_DOUBLE_PRECISION, mixCplRank,
     |       (myAppId+voltId)*100, CplComm, ierr)
          if (ierr /= 0)
     |        call handle_mpi_err(ierr, 'export: TIEGCM send to REMIX')
          endif

          if ( nmixoutapex .ne. 0) then
          call mpi_send(avar2d, nmlat*nmlonp1*nmixoutapex,
     |       MPI_DOUBLE_PRECISION, mixCplRank,
     |       (myAppId+voltId)*100, CplComm, ierr)
          if (ierr /= 0)
     |        call handle_mpi_err(ierr, 'export: TIEGCM send to REMIX')
          endif

        endif

      end subroutine export_remix
!-------------------------------------------------------------------

      subroutine init_mpi_remix(glatp2,glonp1,gmlat,gmlon)

      use params_module,only: nlatp2,nlonp1
      use mpi_module,only: myAppId,voltId,CplComm,mixCplRank,
     |  handle_mpi_err
      use mpi

      integer :: ierr,i
      real(rp),dimension(nlatp2) :: glatp2
      real(rp),dimension(nlonp1) :: glonp1
      real(rp),dimension(nmlat):: gmlat
      real(rp),dimension(nmlonp1):: gmlon

      if (CplComm /= MPI_COMM_NULL) then
        i = 0

        write(*,*) "TGCM: Start init_mpi_remix"
        write(*,*) "TGCM: Sending GEO"
        i = i + 1
        call mpi_send(nlatp2, 1, MPI_INTEGER,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(nlonp1, 1, MPI_INTEGER,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(glatp2, nlatp2, MPI_DOUBLE_PRECISION,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(glonp1, nlonp1, MPI_DOUBLE_PRECISION,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        write(*,*) "TGCM: Sending APEX"
        i = i + 1
        call mpi_send(nmlat, 1, MPI_INTEGER,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(nmlonp1, 1, MPI_INTEGER,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(gmlat, nmlat, MPI_DOUBLE_PRECISION,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        i = i + 1
        call mpi_send(gmlon, nmlonp1, MPI_DOUBLE_PRECISION,
     |    mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        write(*,*) "TGCM: Done init_mpi_remix"
      endif

      end subroutine init_mpi_remix

!-------------------------------------------------------------------

      subroutine import_remix(avar2d,gvar2d)

      use params_module,only: nlatp2,nlonp1
      use mpi_module,only: myAppId,voltId,CplComm,mixCplRank,
     |  handle_mpi_err
      use mpi

        integer :: v,ierr,i,j
        real(rp),dimension(:,:,:) :: avar2d,gvar2d

        if (CplComm /= MPI_COMM_NULL) then
          !write(*,*) "GCMCPL: ", mixCplRank,(myAppId+voltId)*100
          if (nmixingeo .ne. 0) then
            call mpi_recv(gvar2d, nlatp2*nlonp1*nmixingeo,
     |        MPI_DOUBLE_PRECISION, mixCplRank,
     |        (myAppId+voltId)*100, CplComm,
     |        MPI_STATUS_IGNORE, ierr)
            if (ierr /= 0) call handle_mpi_err(ierr,
     |        'import: TIEGCM recv from REMIX')
          endif
          if (nmixinapex .ne. 0) then
            call mpi_recv(avar2d, nmlat*nmlonp1*nmixinapex,
     |        MPI_DOUBLE_PRECISION, mixCplRank,
     |        (myAppId+voltId)*100, CplComm,
     |        MPI_STATUS_IGNORE, ierr)
            if (ierr /= 0) call handle_mpi_err(ierr,
     |        'import: TIEGCM recv from REMIX')
          endif

!!!            call mpi_recv(aurorabc, nmlonp1*2,
!!!     |        MPI_DOUBLE_PRECISION, mixCplRank,
!!!     |        (myAppId+voltId)*100, CplComm,
!!!     |        MPI_STATUS_IGNORE, ierr)
!!!            if (ierr /= 0) call handle_mpi_err(ierr,
!!!     |        'import: TIEGCM recv from REMIX')

        endif

      end subroutine import_remix

!-------------------------------------------------------------------

      subroutine init_mpi_hidra(commId,hidraCplRank)

      use params_module,only: nlat,nlon,nlevp1
      use params_module,only: glon,glat,zpmid
      use cons_module,only  : p0
      use mpi_module,only : CplComm,
     |  handle_mpi_err
      use mpi

      integer :: ierr,icolor,commId,hidraCplRank

      if (CplComm /= MPI_COMM_NULL) then
        write(*,*) "TGCM: Start init_mpi_hidra"

        icolor = 1
        call mpi_send(nlat, 1, MPI_INTEGER,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        icolor = icolor + 1
        call mpi_send(nlon, 1, MPI_INTEGER,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        icolor = icolor + 1
        call mpi_send(nlevp1, 1, MPI_INTEGER,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        ! Convert p0 from "standard pressure" to mbar (hPA)
        icolor = icolor + 1
        call mpi_send(p0/1000., 1, MPI_DOUBLE_PRECISION,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        icolor = icolor + 1
        call mpi_send(glat, nlat, MPI_DOUBLE_PRECISION,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        icolor = icolor + 1
        call mpi_send(glon, nlon, MPI_DOUBLE_PRECISION,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        icolor = icolor + 1
        call mpi_send(zpmid, nlevp1, MPI_DOUBLE_PRECISION,
     |    hidraCplRank, commId*100+icolor, CplComm,ierr)
        if (ierr /= 0)
     |    call handle_mpi_err(ierr, 'initialize: TIEGCM send to REMIX')

        write(*,*) "TGCM: Done init_mpi_hidra"
      endif

      end subroutine init_mpi_hidra
!-----------------------------------------------------------------------
      subroutine finalize_coupling

      if (oneway) then
#ifdef GAMERA
        call ca_close
#endif
      endif

      end subroutine finalize_coupling
!-----------------------------------------------------------------------
      subroutine mage_ucurrent(rim1,rim2,mlon0,mlon1,mlat0,mlat1,nsrhs)
!
! 3/20/14 btf: This routine moved from old serial dynamo (where it was
!   called after sub transf, and followed by 3 calls to mag2geo). 
!   The call to this routine is now in pdynamo after sub complete_integrals, 
!   however as of this date. 
!
! DESCRIPTION:
!
! This subroutine calculate height-integrated neutral wind generated
! field-aligned current (dynamo) to be passed to the M-I coupler to
! solve electric potential. This subroutine is based on the subroutine
! 'nosocoef' in 'current.F' written by Astrid Maute.
! The height-integrated neutral wind field-alined current is calculated in a
! Quisi-Dipole Coordinate that is defined in detail in Richmond (1995). This
! coordinate system removes the 1/|sinI_m| factor in the partial differential
! equation for the electric potential in the Modified Apex Coordinate system.
! 1/|sinI_m| is not defined at magnetic equator in the Modified Apex
! Coordinate system, but is well defined in the Quisi-Dipole Coordinate system
! (I still need to see how this works).
! The neutral dynamo currents are already calculted in subroutine
! 'fieldline-integrals' in the 'dynamo.F' as a global variable
! 'rim(nmlonp1,nmlat,2)'. Subroutine 'rshpde' has the formula to calculate
! height-integrated neutral wind current, but the current there is the sum
! of two hemispheres. We want a global distribution of this current for the M-I
! coupler. Thus the code here is an expanded version of that in "rhspde", but a
! stripped version of "nosocoef". 'nosocoef' also calculates other
! coefficients (lhs) for the potential equation to obtain total field-aligned currents
! including both magnetosphere and thermosphere originated currents. We only need
! thermospheric originated currents for the MAGE M-I coupler.
!
! This subroutine is called by subroutine 'dynamo' after 'call transfer'
! in 'dynamo.F'
!                                                  ---------- Wenbin Wang 09/20/05
! USES
!
      use params_module,only: nmlon,nmlonp1,nmlat,nmlath
      use cons_module,only: dlonm,dlatm,pi,r0
      use mpi_module,only: mp_mag_periodic_f2d

!
! PARAMETERS:
! RETURN VALUE: nsrhs(nmlonp1,nmlat)   ! defined as a global variable above
!
! !REVISION HISTORY:
!
! EOP
!
! Calculate height-integrated field-aligned neutral wind currents for both hemisphere
!
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1
      real(rp),intent(in) :: rim1(mlon0-1:mlon1+1,mlat0-1:mlat1+1),
     |                      rim2(mlon0-1:mlon1+1,mlat0-1:mlat1+1)
      real(rp),intent(out) :: nsrhs(mlon0-1:mlon1+1,mlat0-1:mlat1+1)

!
! Local:
!
      real(rp),parameter :: unitvm(nmlon)=1.
      real(rp) :: cs(mlat0-1:mlat1+1)
      real(rp) :: dfac
      real(rp) :: fmsub(mlon0:mlon1,mlat0:mlat1)
      integer :: j,i
      integer mlon00,mlon11,mlat00,mlat11
      integer latstart, latend
!
! Externals:
!
      real,external :: sddot ! in util.F
!
! Set halos indexes
!
      mlon00=mlon0-1 ; mlon11=mlon1+1
      mlat00=mlat0-1 ; mlat11=mlat1+1
!
! Calculate coefficients for dynamo pde for both hemisphere
!
!
! Clear arrays
!
      nsrhs(:,:)    = 0.0
!
! Calculate magnetic latitude cosin array
!
      do j = mlat00,mlat11     
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!
! Calculate right hand side of pde from rim1 and rim2
!
      latstart=mlat0
      latend=mlat1
      if (latstart==1) then 
        latstart=2
      endif
      if (latend==nmlat) then 
        latend=nmlat-1
      endif
      do j = latstart,latend  
        if (j == nmlath) cycle  ! skip equator
        !
        ! Differentiate rim(1) w.r.t lamda
        !
        do i = mlon0,mlon1
          if (i /= 1 .and. i /= nmlonp1) then  ! Skip longitudinal boundaries
            nsrhs(i,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim1(i+1,j)-rim1(i-1,j))
          else ! Values at longitudinal boundaries
            if (i==1) then
              nsrhs(1,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim1(2,j)-rim1(0,j))
            else
              nsrhs(nmlonp1,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim1(mlon0,j)-rim1(nmlonp1-1,j))
            endif
          endif
        enddo
      enddo  ! j loop
!
! Differentiate rim(2) w.r.t theta0
!
      do j = mlat0,mlat1  ! south pole to north pole
        if (j==1) then ! Calculate south pole value
           nsrhs(mlon0,1)     = -2./float(nmlon)*
     |               sddot(nmlon,unitvm,rim2(mlon0,2))/cs(2)
           ! Extend over longitude
           nsrhs(:,1) = nsrhs(mlon0,1)
        elseif (j==nmlat) then ! Calculate calculate north pole value 
           nsrhs(mlon0,nmlat) = -2./float(nmlon)*
     |               sddot(nmlon,unitvm,rim2(mlon0,nmlat-1))/cs(mlat0-1)
           ! Extend over longitude
           nsrhs(:,nmlat) = nsrhs(mlon0,nmlat)
         elseif (j==nmlath) then ! Calculate equator values
           do i=mlon0,mlon1
             if (i==1) then
               nsrhs(i,j) = 0.5/dlonm*(rim1(i+1,j)-rim1(i-1,j))
               nsrhs(i,j) = nsrhs(i,j) + 1./dlatm*(cs(j)*
     |                 rim2(i,j)+ cs(j+1)*rim2(i,j+1))
             elseif (i==nmlonp1) then
               nsrhs(i,j) = 0.5/dlonm*(rim1(nmlonp1+1,j)-rim1(i-1,j))
               nsrhs(i,j) = nsrhs(i,j) + 1./dlatm*(cs(j)*
     |                 rim2(i,j)+ cs(j+1)*rim2(i,j+1))
             else
               nsrhs(i,j) = 0.5/dlonm*(rim1(i+1,j)-rim1(i-1,j))
               nsrhs(i,j) = nsrhs(i,j) + 1./dlatm*(cs(j)*
     |                    rim2(i,j)+ cs(j+1)*rim2(i,j+1))
             endif

           enddo

         elseif (j<nmlath) then !calculate south pole to equator
           do i = mlon0,mlon1
             nsrhs(i,j)  = nsrhs(i,j) - 1.0/(dlatm*cs(j))*0.5*
     |               (rim2(i,j+1)*cs(j+1)-rim2(i,j-1)*cs(j-1))
           enddo
         else !calculate north pole to equator
           do i = mlon0,mlon1
             nsrhs(i,j) = nsrhs(i,j) + 1.0/(dlatm*cs(j))*0.5*
     |               (rim2(i,j+1)*cs(j+1)-rim2(i,j-1)*cs(j-1))
           enddo
         endif
      enddo
!
! Periodic points
!
      fmsub(:,:)=nsrhs(mlon0:mlon1,mlat0:mlat1)
      call mp_mag_periodic_f2d(fmsub,mlon0,mlon1,mlat0,mlat1,1)

      nsrhs(mlon0:mlon1,mlat0:mlat1) = fmsub(:,:)
!
! Scale rhs by refernce radius (R_E + H0) in meters dfac = r0*1e-2
!
      dfac = r0*1.0e-2
      nsrhs(:,:) = -1.*nsrhs(:,:)/dfac
!
      end subroutine mage_ucurrent


      end module mage_coupling_module
