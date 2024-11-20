
#ifdef GAMERA
      module mage_oneway
      use input_module,only: mixfile,mxlen_filename
      use params_module,only: nlat, nlon, rp
      use nchist_module,only:handle_ncerr
      use dates
      use netcdf

      real*8, allocatable,dimension(:,:,:,:) :: mix_interpolant
      real*8, allocatable,dimension(:,:,:,:) :: tgcm_interpolant
      real*8, allocatable,dimension(:,:) :: tlat, tlon,
     |  tcolat,tlon2,tlat2
      real*8, allocatable,dimension(:,:) :: xtgcm,ytgcm,ztgcm
      real*8, allocatable,dimension(:,:) :: mixX,mixY,mixZ
      real*8, allocatable,dimension(:,:) :: mixlat,mixlon,mixcolat
      real*8, allocatable,dimension(:,:) :: mixPotN,mixEngN,
     |  mixFluxN,mixBEngN
      real*8, allocatable,dimension(:,:) :: mixBFluxN,mixCEngN,mixCFluxN
      real*8, allocatable,dimension(:,:) :: mixPotS,mixEngS,
     |  mixFluxS,mixBEngS
      real*8, allocatable,dimension(:,:) :: mixBFluxS,mixCEngS,mixCFluxS
      real*8, allocatable,dimension(:,:) :: tmpvar
      real*8, allocatable,dimension(:,:) :: mixOutX,mixOutY,mixOutZ
      real*8 :: LowLatBoundary
      real*8, parameter :: TINY=1.e-8
      real*8,dimension(:,:,:),allocatable :: ting_ped_interp,
     |  ting_hall_interp
      real(rp),allocatable,dimension(:,:) :: mixaurorabc1d
      integer :: nmixlon,nmixlat
      integer,parameter :: iSMtoGSM=1,iGSMtoSM=-1,
     |  iGEOtoGSM=1,iGSMtoGEO=-1
      real*8, allocatable, dimension(:) :: MJDs
      integer :: s0, sE, nSteps,cStep

      real*8, dimension(:,:),allocatable ::
     |  ting_pot_interp,  ! potential in geo coord from MIX
     |  ting_eng_interp,  ! energy in geo coord from MIX
     |  ting_flux_interp, ! flux in geo coord from MIX
     |  ting_bbeAvgEng_interp, ! BBE Average Energy geo coord from MIX
     |  ting_bbeNumFlux_interp,  ! BBE Number flux geo coord from MIX
     |  ting_cuspAvgEng_interp, ! Cusp Average Energy geo coord from MIX
     |  ting_cuspNumFlux_interp  ! Cusp Number flux geo coord from MIX

      integer :: h5fId,inId,gId

      integer :: nlatp2,cplStep = 1
      logical :: do_init = .true., firstread = .true.
      logical :: step1140Fix = .true.
      real*8, parameter :: pi = 4*atan(1.0d0)


      contains
!-----------------------------------------------------------------------
!
! Functions for using ad hoc coupling in CISM_LTR
!
!-------------------------------------------------------------------
      subroutine ca_init
!
! This subroutine initialize ad hoc interface for data exchange
!
! Open Fortran I/O units for M-I coupling outputs, will be replaced with
! netcdf calls later
!
      use params_module,only: glon,glat,nlat,nlonp1
      use mpi_module,only: mytid,TIEGCM_WORLD
      use mpi
!
! Local:
!
      real*8,dimension(4,4) :: A,Ainv
      real*8 :: dp,tempLat,tempLon
      real*8,allocatable,dimension(:,:) :: gridX,gridY
      integer :: i,j,jp1,nji(2)
      character(len=mxlen_filename) :: filename
      !character(len=240) :: mixfile = "msphere.mix.h5"
      integer :: npath
      logical isFileOK


! Get MIX grid information to initialize reader variables
      !call getfile(mix_restartfile,filename)
! read grid

      ! Allocate things
      if (.not. allocated(ting_pot_interp)) then
        allocate(ting_pot_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_eng_interp)) then
        allocate(ting_eng_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_flux_interp)) then
        allocate(ting_flux_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_bbeAvgEng_interp)) then
        allocate(ting_bbeAvgEng_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_bbeNumFlux_interp)) then
        allocate(ting_bbeNumFlux_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_cuspAvgEng_interp)) then
        allocate(ting_cuspAvgEng_interp(nlat,nlonp1))
      endif
      if (.not. allocated(ting_cuspNumFlux_interp)) then
        allocate(ting_cuspNumFlux_interp(nlat,nlonp1))
      endif

      if (mytid == 0) then
      filename = ' '
      call getfile(mixfile,filename)
      write(6,*) 'inside ca_init:',trim(filename)
      call openMixFile(trim(filename))
      !call CheckAndKill(trim(mix2gcmLock))
      write(6,*) "TGCM: Stepinfo",mytid
      call StepInfo(trim(filename),s0,sE,nSteps)
      write(6,*) "TGCM: nSteps",mytid,s0,sE,nSteps
      if (.not. allocated(MJDs) ) allocate(MJDs(nSteps))
      call StepMJDs(trim(filename),s0,sE,MJDs)
      write(6,*) "TGCM: StepMJDs",mytid,MJDs
      write(6,*) "TGCM init done",mytid
      cStep = 1

      !call readDimension(ncid,'phony_dim_0',nmixlat)
      !call readDimension(ncid,'phony_dim_1',nmixlon)

      !if (.not.allocated(gridX)) then
      !    allocate(gridX(nmixlat,nmixlon))
      !end if
      !
      !if (.not.allocated(gridY)) then
      !    allocate(gridY(nmixlat,nmixlon))
      !end if
      !
      call closeGroup()
      call read2DVariable("X",gridX)
      call read2DVariable("Y",gridY)

      !size of the Mix grid
      nji=shape(gridX)
      nmixlat = nji(1)-1 !t
      nmixlon = nji(2)-1 !p
      write(6,*) "read grid size:", nmixlat,nmixlon

      endif

      call mpi_bcast(s0,1,MPI_INTEGER,0,
     |  TIEGCM_WORLD,ierr)
      call mpi_bcast(sE,1,MPI_INTEGER,0,
     |  TIEGCM_WORLD,ierr)
      call mpi_bcast(nSteps,1,MPI_INTEGER,0,
     |  TIEGCM_WORLD,ierr)
      call mpi_bcast(nmixlat,1,MPI_INTEGER,0,
     |  TIEGCM_WORLD,ierr)
      call mpi_bcast(nmixlon,1,MPI_INTEGER,0,
     |  TIEGCM_WORLD,ierr)
      if (.not. allocated(gridX))
     |  allocate(gridX(nmixlat+1,nmixlon+1))
      if (.not. allocated(gridY))
     |  allocate(gridY(nmixlat+1,nmixlon+1))
      call mpi_bcast(gridX,(nmixlat+1)*(nmixlon+1),MPI_REAL8,0,
     |  TIEGCM_WORLD,ierr)
      call mpi_bcast(gridY,(nmixlat+1)*(nmixlon+1),MPI_REAL8,0,
     |  TIEGCM_WORLD,ierr)



! Allocatable MIX variable arrays
      if (.not.allocated(mixX)) then
         allocate(mixX(nmixlon,nmixlat))
      endif
      if (.not.allocated(mixY)) then
         allocate(mixY(nmixlon,nmixlat))
      endif
      if (.not.allocated(mixZ)) then
         allocate(mixZ(nmixlon,nmixlat))
      endif
      if (.not.allocated(mixlat)) then
         allocate(mixlat(nmixlon,nmixlat))
      endif
      if (.not.allocated(mixlon)) then
         allocate(mixlon(nmixlon,nmixlat))
      endif
      ! Allocate North Arrays
      if (.not.allocated(mixPotN)) then
         allocate(mixPotN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixEngN)) then
         allocate(mixEngN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixFluxN)) then
         allocate(mixFluxN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixBEngN)) then
         allocate(mixBEngN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixBFluxN)) then
         allocate(mixBFluxN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixCEngN)) then
         allocate(mixCEngN(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixCFluxN)) then
         allocate(mixCFluxN(nmixlat,nmixlon))
      endif
      ! Allocate South Arrays
      if (.not.allocated(mixPotS)) then
         allocate(mixPotS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixEngS)) then
         allocate(mixEngS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixFluxS)) then
         allocate(mixFluxS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixBEngS)) then
         allocate(mixBEngS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixBFluxS)) then
         allocate(mixBFluxS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixCEngS)) then
         allocate(mixCEngS(nmixlat,nmixlon))
      endif
      if (.not.allocated(mixCFluxS)) then
         allocate(mixCFluxS(nmixlat,nmixlon))
      endif

      if (.not. allocated(mixaurorabc1d)) then
         allocate(mixaurorabc1d(nmixlon,2))
      endif

! Initialize MIX variables
      do i=1,nmixlat
         do j=1,nmixlon
            mixX(j,i) = gridX(i,j)
            mixY(j,i) = gridY(i,j)
            mixZ(j,i) = sqrt(1. -gridX(i,j)*gridX(i,j)-
     |  gridY(i,j)*gridY(i,j))
            mixlat(j,i) = asin(mixZ(j,i))
            mixlon(j,i) = atan2(mixY(j,i),mixX(j,i))
            if (mixlon(j,i).lt.0.0) mixlon(j,i)=mixlon(j,i)+2*pi
         enddo
      enddo
      mixlon(:,1) = mixlon(:,2)
      LowLatBoundary = minval(mixlat)
      ! North
      mixPotN   = 0.0
      mixEngN   = 0.0
      mixFluxN  = 0.0
      mixBEngN  = 0.0
      mixBFluxN = 0.0
      mixCEngN  = 0.0
      mixCFluxN = 0.0
      ! South
      mixPotS   = 0.0
      mixEngS   = 0.0
      mixFluxS  = 0.0
      mixBEngS  = 0.0
      mixBFluxS = 0.0
      mixCEngS  = 0.0
      mixCFluxS = 0.0
      if (mytid == 0) then
         !write(6,*) "mix Z: ",mixZ
         write(6,*) "mix Lat grid1: ", mixlat(1,:)
         write(6,*) "mix Lat grid2: ", mixlat(nmixlon-1,:)
         write(6,*) "mix Lat grid3: ", mixlat(nmixlon,:)
         write(6,*) "mix Lon grid1: ", mixlon(:,1)
         write(6,*) "mix Lon grid2: ", mixlon(:,2)
         write(6,*) "mix LowLatBoundary: ",LowLatBoundary
         write(6,*) "mix pi: ",pi
      endif

      if (.not. allocated(mixOutX)) allocate(mixOutX(nmixlat+1,nmixlon))
      if (.not. allocated(mixOutY)) allocate(mixOutY(nmixlat+1,nmixlon))
      if (.not. allocated(mixOutZ)) allocate(mixOutZ(nmixlat+1,nmixlon))

      write(6,*) " Done Allocating"
      !Figure out which end of the latitude range we are attaching to.
      if (mixlat(1,1) < 7./8.*pi/2.) then
        do j = 1,nmixlon
          tempLat = mixlat(j,1) - (mixlat(j,2) - mixlat(j,1))
          mixOutZ(1,j) = sin(tempLat)
          mixOutY(1,j) = cos(tempLat)*sin(mixlon(j,1))
          mixOutX(1,j) = cos(tempLat)*cos(mixlon(j,1))
          do i = 1,nmixlat
            mixOutZ(i+1,j) = mixZ(j,i)
            mixOutY(i+1,j) = mixY(j,i)
            mixOutX(i+1,j) = mixX(j,i)
          end do
        end do
      else
        do j = 1,nmixlon
          tempLat = mixlat(j,nmixlat) - (mixlat(j,nmixlat-1) -
     |  mixlat(j,nmixlat))
          mixOutZ(nmixlat+1,j) = sin(tempLat)
          mixOutY(nmixlat+1,j) = cos(tempLat)*sin(mixlon(j,1))
          mixOutX(nmixlat+1,j) = cos(tempLat)*cos(mixlon(j,1))
          do i = 1,nmixlat
            mixOutZ(i,j) = mixZ(j,i)
            mixOutY(i,j) = mixY(j,i)
            mixOutX(i,j) = mixX(j,i)
          end do
        end do
       end if

      write(6,*) " Done expanding MixOut"
! Initialize MIX interpolant
      if (.not.allocated(mix_interpolant)) then
         allocate(mix_interpolant(nmixlon,nmixlat,4,4))
      endif
      write(6,*) " Start Interpolant"
      mix_interpolant = 0.0
! Construct MIX interpolant
      do i=2,nmixlat !pole?
         do j=1,nmixlon !pole?
            if (j.ne.nmixlon) then
               jp1=j+1
               dp = 0.0
            else
               jp1=1
               dp = 2*pi
            endif

            A = reshape(
     |       (/
     |       1.0,mixlon(j,i),mixlat(j,i),mixlon(j,i)*mixlat(j,i),
     |        1.0,mixlon(jp1,i)+dp,mixlat(jp1,i),
     |            (mixlon(jp1,i)+dp)*mixlat(jp1,i),
     |      1.0,mixlon(j,i-1),mixlat(j,i-1),mixlon(j,i-1)*mixlat(j,i-1),
     |        1.0,mixlon(jp1,i-1)+dp,mixlat(jp1,i-1),
     |            (mixlon(jp1,i-1)+dp)*mixlat(jp1,i-1)
     |        /)
     |        ,shape(A))

            ! exclude the poles to avoid inverting the bad matrix there

            if (.not.(
     |         ( (i.eq.1).and.(pi/2-mixlat(1,1).le.1.0e-8) )
     |         )) then
               call matinv4(A,Ainv)
               mix_interpolant(j,i,:,:) = Ainv ! for cell (j,i) or (Lon,Lat)
            endif
         enddo
      enddo

      write(6,*) " Done Interpolant"
! Initialize the global tgcm search array
      nlatp2 = nlat+2
      if (.not. allocated(tlat)) allocate(tlat(nlon,nlat),
     |  tcolat(nlon,nlat),tlon(nlon,nlat))
      if (.not. allocated(tlat2)) allocate(tlat2(nlon,nlatp2),
     |  tlon2(nlon,nlatp2))
      if (.not. allocated(xtgcm)) allocate(xtgcm(nlon,nlat),
     |  ytgcm(nlon,nlat),ztgcm(nlon,nlat))
      if (.not. allocated(tgcm_interpolant)) allocate(
     |  tgcm_interpolant(nlon,nlatp2,4,4))
      write(6,*) "ALLOCATE TGCM GRID"
      do i = 1,nlat
         do j = 1,nlon
            tlat(j,i) = glat(i)*pi/180.
            tlon(j,i) = glon(j)*pi/180.
            tcolat(j,i) = (90.-glat(i))*pi/180.
            if (tlon(j,i) .lt. 0.0) tlon(j,i) = tlon(j,i) + 2*pi
            xtgcm(j,i) = cos(tlat(j,i))*cos(tlon(j,i))
            ytgcm(j,i) = cos(tlat(j,i))*sin(tlon(j,i))
            ztgcm(j,i) = sin(tlat(j,i))
            tlat2(j,i+1) = glat(i)*pi/180.
            tlon2(j,i+1) = glon(j)*pi/180.
         enddo
      enddo
      write(*,*) "TGCM GRID EXPAND1"
      tlat2(:,1) = -pi/2.
      write(*,*) "TGCM GRID EXPAND2"
      tlat2(:,nlatp2) = pi/2.
      write(*,*) "TGCM GRID EXPAND3"
      tlon2(:,1) = tlon2(:,2)
      write(*,*) "TGCM GRID EXPAND4"
      tlon2(:,nlatp2) = tlon2(:,nlat)
      if (mytid == 0) then
      write(6,*) "glatt: ",glat
      write(6,*) "glont: ",glon
      write(6,*) "tcolat: ",tcolat(1,:)
      write(6,*) "tlat: ", tlat(1,:)
      write(6,*) "tlon: ", tlon(:,2)
      endif

! Initialize TGCM interpolant
      tgcm_interpolant = 0.
      ! note, looping over cells -- not nodes -- here
      do i=2,nlat+1
         do j=1,nlon  ! we go to +1 point in phi to cover the gap between Np and 1
            if (j.ne.nlon) then
               jp1=j+1
               dp = 0.
            else
               jp1=1
               dp = 2*pi
            endif

            A = reshape(
     |            (/
     |            1.0,tlon2(j,i),       tlat2(j,i),
     |               tlon2(j,i)        *tlat2(j,i),
     |            1.0,tlon2(jp1,i)  +dp,tlat2(jp1,i),
     |               (tlon2(jp1,i)  +dp)*tlat2(jp1,i),
     |            1.0,tlon2(j,i+1),     tlat2(j,i+1),
     |               tlon2(j,i+1)      *tlat2(j,i+1),
     |            1.0,tlon2(jp1,i+1)+dp,tlat2(jp1,i+1),
     |               (tlon2(jp1,i+1)+dp)*tlat2(jp1,i+1)
     |            /)
     |            ,shape(A))

            ! exclude the poles to avoid inverting the bad matrix there
            if (.not.(
     |            ( (i.eq.1).and.(pi/2-abs(
     |               tlat2(1,1)).le.TINY) ) .or.
     |            ( (i.eq.nlat+2).and.((pi/2-abs(
     |               tlat2(1,nlatp2))).le.TINY) )
     |            )) then
               call matinv4(A,Ainv)
               tgcm_interpolant(j,i,:,:) = Ainv   ! for cell (j,i)
            end if
         enddo
      enddo

      if (.not. allocated(ting_ped_interp))
     |   allocate(ting_ped_interp(nmixlat+1,nmixlon,2))
      if (.not. allocated(ting_hall_interp))
     |   allocate(ting_hall_interp(nmixlat+1,nmixlon,2))

!
!     End initialize ad hoc coupling
!
      write(6,*) 'Ad hoc coupling Init Complete'

      end subroutine ca_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_receive(gpot,geng,gflx,gpotm,aurllbN,aurllbS)
!
! This is a subroutine to receive auroral data via file exchanges
! via HDF5 format from REMIX
!
      use params_module,only: nlat,nlon,glat,nlonp4,nlatp1,nlonp1
      use init_module,only: iyear,iday,secs
      use mpi_module,only: mytid,TIEGCM_WORLD
      use nchist_module,only:handle_ncerr
      use cons_module,only: rtd
      use geopack
      use mpi

      real,dimension(nlonp4,nlat) ::
     |  gpot,        ! potential in geographic coordinates, periodic boundary
     |  geng,        ! energy in geographic coordinates, periodic boundary
     |  gflx         ! flux in geographic coordinates, periodic boundary
      real,dimension(nlonp4,0:nlatp1) :: gpotm

!
! Local:
!
      character(len=mxlen_filename) :: filename
      integer,dimension(5) :: utime
      real*8,dimension(4) :: x,mapping
      real*8,dimension(4) :: FPot,FEng,FFlux,FBEng,FBFlux,FCEng,FCFlux
      real*8 :: xgsm,ygsm,zgsm,xsm,ysm,zsm
      real*8 :: tempLat, tempLon,tempproduct
      real*8 :: tempcolat
      real*8 :: mixtime,gcmtime,csec
      integer :: cyr,cdoy,cmon,cday,chr,cmin
      integer :: i,j,j1,i1,j1p1
      integer :: istat, oldStep, req(7)
      integer :: id_nlat,id_nlon,idv_colat
      integer :: idv_lon,idv_pot,idv_flx,idv_eng
      character(len=240) :: char240,gStr
      logical :: fExist = .false.
      real(rp) :: nfluxllb = 1.0e7_rp

      if (cplStep == 1 .and. do_init) then
        call ca_init
        do_init = .false.
      end if

      if (mytid == 0) then
      write(6,*) "ca_receive",firstread

      !inquire(file=mix2gcmLock,exist=fExist)
      !
      !write(*,*) "TGCM: waiting for ",trim(mix2gcmLock)," to appear"
      !do while (.not.fExist)
      !   inquire(file=mix2gcmLock,exist=fExist)
      !   call sleep(1)
      !end do
      !write(*,*) "TGCM: ",trim(mix2gcmLock)," is gone."
      !write(6,*) "TGCM: reading ",trim(mix2gcmH5),' start coupling @ ',cplStep
      !call CheckAndKill(cplStr)

      filename = ' '
      call getfile(mixfile,filename)

      utime(1) = iyear
      utime(2) = iday
      utime(3) = floor(secs/3600.)
      utime(4) = floor((secs - utime(3)*3600.)/60.)
      utime(5) = secs - utime(3)*3600. - utime(4)*60.


      call mjd2ut(MJDs(cStep),cyr,cdoy,cmon,cday,chr,cmin,csec)
      !call ut2mjd(cMJD,utime(1),0,utime(2),utime(3)

      oldStep = cStep

      gcmtime = utime(2) + secs/(24.*60.*60.)
      mixtime = cdoy+chr/24.+cmin/(24.*60.)+csec/(24.*60.*60.)
      do while(mixtime < gcmtime)
        if (cStep == nSteps) then
            exit
        endif
        cStep = cStep + 1
        call mjd2ut(MJDs(cStep),cyr,cdoy,cmon,cday,chr,cmin,csec)
        mixtime = cdoy+chr/24.+cmin/(24.*60.)+csec/(24.*60.*60.)
      end do

      if (mixtime > gcmtime) then
        cStep = max(cStep-1,1)
      endif


      if (oldStep .ne. cStep) then
        firstread = .true.
      endif

      write(*,*) "MIXREAD:",mixtime,gcmtime,cStep,oldStep,firstread

      if (firstread) then
      firstread = .false.

      ting_pot_interp = 0.0
      ting_eng_interp = 0.0
      ting_flux_interp = 0.0
      ting_bbeAvgEng_interp = 0.0
      ting_bbeNumFlux_interp = 0.0
      ting_cuspAvgEng_interp = 0.0
      ting_cuspNumFlux_interp = 0.0

      write(gStr,'(A,I0)') "Step#", s0 + cStep - 1
      write(6,*) "REMIX READING: ", cStep, MJDs(cStep),cyr,
     |  cdoy,cmon,cday,chr,cmin,csec

      !call openMixFile(filename,gStr)
      call openGroup(gStr)
      call read2DVariable("Potential NORTH",mixPotN)
      call read2DVariable("Number flux NORTH",mixFluxN)
      call read2DVariable("Average energy NORTH",mixEngN)

      ! SOUTH
      call read2DVariable("Potential SOUTH",mixPotS)
      call read2DVariable("Number flux SOUTH",mixFluxS)
      call read2DVariable("Average energy SOUTH",mixEngS)

      if (mytid==0) then
      write(6,*) "TGCM: before transformN: ", minval(mixPotN),
     |  maxval(mixPotN)
      write(6,*) "TGCM: before transformS: ", minval(mixPotS),
     |  maxval(mixPotS)
      endif

      mixPotN = mixPotN*1.e3 !Volts!
      mixBEngN  = 0.
      mixBFluxN = 0.
      mixCEngN  = 0.
      mixCFluxN = 0.
      mixPotS = mixPotS*1.e3 !Volts!
      mixBEngS  = 0.
      mixBFluxS = 0.
      mixCEngS  = 0.
      mixCFluxS = 0.

      ! find auroral boundary
      do j=1,nmixlon
         do i=1,nmixlat
            if (mixfluxN(i,j) >= nfluxllb) exit
         enddo
         i = min(i,nmixlat)
         mixaurorabc1d(j,1) = max(90.0_rp-mixlat(j,i)*rtd,15._rp) ! aurllbj is co-lat.
      enddo
      ! find auroral boundary
      do j=1,nmixlon
         do i=1,nmixlat
            if (mixfluxS(i,j) >= nfluxllb) exit
         enddo
         i = min(i,nmixlat)
         mixaurorabc1d(j,2) = max(90.0_rp-mixlat(j,i)*rtd,15._rp) ! aurllbj is co-lat.
      enddo

      endif !firstread

! NORTH
      write(6,*) "TGCM: CA_receive time: ",utime(1),
     | utime(2),utime(3),utime(4),utime(5)
      call RECALC(utime(1),utime(2),utime(3),utime(4),utime(5))
      if (mytid == 0) then
      write(6,*) "TGCM: before PotN: ", minval(mixPotN),maxval(mixPotN)
      write(6,*) "TGCM: minloc: ", minloc(abs(mixPotN-
     |  minval(mixPotN)),2)
      write(6,*) "TGCM: before EngN: ", minval(mixEngN),
     |  maxval(mixEngN)
      write(6,*) "TGCM: before FluxN: ", minval(mixFluxN),
     |  maxval(mixFluxN)
      write(6,*) "TGCM: before BEngN: ", minval(mixBEngN),
     |  maxval(mixBEngN)
      write(6,*) "TGCM: before BFluxN: ",
     |  minval(mixBFluxN),maxval(mixBFluxN)
      write(6,*) "TGCM: before CEngN: ", minval(mixCEngN),
     |  maxval(mixCEngN)
      write(6,*) "TGCM: before CFluxN: ",
     |  minval(mixCFluxN),maxval(mixCFluxN)
      write(6,*) "TGCM: shape mixPotN:", shape(mixPotN)
      endif

      do i=1,nlat
         do j=1,nlon
            !if (glat(i) .gt. 0) then

      call GEOGSW_08(xtgcm(j,i),ytgcm(j,i),ztgcm(j,i),
     |  xgsm,ygsm,zgsm,iGEOtoGSM)
      call SMGSW_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,iGSMtoSM)
      tempLat = asin(zsm) ! latitude (-pi/2, pi/2)
      tempLon = atan2(ysm,xsm)
      !write(6,*) "tgrid: ",templon,templat
      if (tempLon.lt.0.0) tempLon=tempLon+2*pi
      if (tempLat.le.LowLatBoundary) then
         ! you're below the MIX boundary. Zero out for now
         cycle
      endif
      ! Search on MIX grid for point nearest TGCM rotated point
      call tgcm_search(mixLon,mixLat,tempLon,tempLat,j1,i1)
      if (mytid == 0) then
      !write(6,*) "search: ", tempLon,tempLat,mixLon(j1,i1),mixLat(j1,i1),j1,i1
      !write(6,*) "grid search2: ", j1,i1
      endif
      ! One point mapping
      x = (/1.0,tempLon,tempLat,tempLat*tempLon /)
      if ((i1.eq.1)) then!.and. (pi/2-mixlat(1,1)).lt.1.0e-8) then
         mapping = 1./3.*(/0.5,0.5,1.,1./)
         write(6,*) "TGCM: grid search: ", tempLon,tempLat,
     |     mixLon(j1,i1),mixLat(j1,i1)
         write(6,*) "TGCM: grid search2: ", j1,i1
      else
         mapping = matmul(mix_interpolant(j1,i1,:,:),x)
      endif

      ! check periodicity
      if (j1.eq.nmixlon) then
         j1p1=1
      else
         j1p1=j1+1
      endif

      ! interpolate MIX data to temp (tgcm) position
      FPot = (/mixPotN(i1,j1),mixPotN(i1,j1p1),
     |    mixPotN(i1-1,j1),mixPotN(i1-1,j1p1) /)
      FEng = (/mixEngN(i1,j1),mixEngN(i1,j1p1),
     |    mixEngN(i1-1,j1),mixEngN(i1-1,j1p1) /)
      FFlux = (/mixFluxN(i1,j1),mixFluxN(i1,j1p1),
     |    mixFluxN(i1-1,j1),mixFluxN(i1-1,j1p1) /)
      FBEng = (/mixBEngN(i1,j1),mixBEngN(i1,j1p1),
     |    mixBEngN(i1-1,j1),mixBEngN(i1-1,j1p1) /)
      FBFlux = (/mixBFluxN(i1,j1),mixBFluxN(i1,j1p1),
     |    mixBFluxN(i1-1,j1),mixBFluxN(i1-1,j1p1) /)
      FCEng = (/mixCEngN(i1,j1),mixCEngN(i1,j1p1),
     |    mixCEngN(i1-1,j1),mixCEngN(i1-1,j1p1) /)
      FCFlux = (/mixCFluxN(i1,j1),mixCFluxN(i1,j1p1),
     |    mixCFluxN(i1-1,j1),mixCFluxN(i1-1,j1p1) /)

      ting_pot_interp(i,j) = dot_product(mapping,FPot)
      tempproduct = mapping(1)*FPot(1)+mapping(2)*FPot(2)+
     |    mapping(3)*Fpot(3)+mapping(4)*FPot(4)
      ting_eng_interp(i,j) = max(dot_product(mapping,FEng),0.0)
      ting_flux_interp(i,j) = max(dot_product(mapping,FFlux),0.0)
      ting_bbeAvgEng_interp(i,j) = max(dot_product(mapping,FBEng),0.0)
      ting_bbeNumFlux_interp(i,j) = max(dot_product(mapping,FBFlux),0.0)
      ting_cuspAvgEng_interp(i,j) = max(dot_product(mapping,FCEng),0.0)
      ting_cuspNumFlux_interp(i,j) = max(dot_product(
     |   mapping,FCFlux),0.0)
      if (mytid == 0 .and. tempproduct.ne.ting_pot_interp(i,j)) then
         write(6,*) "TGCM: OOB POTN: ", ting_pot_interp(i,j),i,j,i1,j1
         write(6,*) "TGCM: OOB LOCN: ", tempLon,tempLat,mixLon(j1,i1),
     |    mixLat(j1,i1)
         write(6,*) "TGCM: OOB VAL: ", FPOT, " mapping:", mapping
      endif
      if (mytid == 0 .and. ting_flux_interp(i,j) .lt. 0.0) then
         write(6,*) "TGCM: OOB Flux: ", ting_flux_interp(i,j),i,j,i1,j1
         write(6,*) "TGCM: OOB LOCN: ", tempLon,tempLat,mixLon(j1,i1),
     |    mixLat(j1,i1)
         write(6,*) "TGCM: OOB VAL: ", FEng
      endif
            !endif
         enddo
      enddo
      if (mytid == 0) then
      write(6,*) "TGCM: mix: min/max potN: ",minval(ting_pot_interp),
     |  maxval(ting_pot_interp)
      write(6,*) "TGCM: mix: min/max engN: ",minval(ting_eng_interp),
     |  maxval(ting_eng_interp)
      write(6,*) "TGCM: mix: min/max fluxN: ",minval(ting_flux_interp),
     |  maxval(ting_flux_interp)
      write(6,*) "TGCM: mix: min/max bengN: ",
     |  minval(ting_bbeAvgEng_interp),
     |  maxval(ting_bbeAvgEng_interp)
      write(6,*) "TGCM: mix: min/max bfluxN: ",
     |  minval(ting_bbeNumFlux_interp),
     |  maxval(ting_bbeNumFlux_interp)
      write(6,*) "TGCM: mix: min/max cengN: ",
     |  minval(ting_cuspAvgEng_interp),
     |  maxval(ting_cuspAvgEng_interp)
      write(6,*) "TGCM: mix: min/max cfluxN: ",
     |  minval(ting_cuspNumFlux_interp),
     |  maxval(ting_cuspNumFlux_interp)
      endif

! SOUTH
      do i=1,nlat
         do j=1,nlon
            !if (glat(i) .lt. 0) then

      call GEOGSW_08(xtgcm(j,i),ytgcm(j,i),ztgcm(j,i),
     |  xgsm,ygsm,zgsm,iGEOtoGSM)
      call SMGSW_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,iGSMtoSM)
      tempLat = asin(zsm) *-1.! latitude (-pi/2, pi/2)
      tempLon = atan2(ysm*-1.,xsm)
      if (tempLon.lt.0.0) tempLon=tempLon+2*pi
      if (tempLat.le.(LowLatBoundary)) then
         ! you're below the MIX boundary. Zero out for now
         cycle
      endif
      ! Search on MIX grid for point nearest TGCM rotated point
      call tgcm_search(mixLon,mixLat,tempLon,tempLat,j1,i1)
      !write(6,*) "grid search3: ", tempLon,tempLat,mixLon(j1,i1),mixLat(j1,i1)*-1.
      !write(6,*) "grid search4: ", j1,i1
      ! One point mapping
      x = (/1.0d0,tempLon,tempLat,tempLat*tempLon /)
      if ((i1.eq.1))then !.and. (pi/2-mixlat(1,1)).lt.1.0e-8) then
         mapping = 1./3.*(/0.5,0.5,1.,1./)
      else
         mapping = matmul(mix_interpolant(j1,i1,:,:),x)
      endif

      ! check periodicity
      if (j1.eq.nmixlon) then
         j1p1=1
      else
         j1p1=j1+1
      endif

      !mapping = 1./4.*(/1.,1.,1.,1./)
      ! interpolate MIX data to temp (tgcm) position
      FPot = (/mixPotS(i1,j1),mixPotS(i1,j1p1),
     |    mixPotS(i1-1,j1),mixPotS(i1-1,j1p1) /)
      FEng = (/mixEngS(i1,j1),mixEngS(i1,j1p1),
     |    mixEngS(i1-1,j1),mixEngS(i1-1,j1p1) /)
      FFlux = (/mixFluxS(i1,j1),mixFluxS(i1,j1p1),
     |    mixFluxS(i1-1,j1),mixFluxS(i1-1,j1p1) /)
      FBEng = (/mixBEngS(i1,j1),mixBEngS(i1,j1p1),
     |    mixBEngS(i1-1,j1),mixBEngS(i1-1,j1p1) /)
      FBFlux = (/mixBFluxS(i1,j1),mixBFluxS(i1,j1p1),
     |    mixBFluxS(i1-1,j1),mixBFluxS(i1-1,j1p1) /)
      FCEng = (/mixCEngS(i1,j1),mixCEngS(i1,j1p1),
     |    mixCEngS(i1-1,j1),mixCEngS(i1-1,j1p1) /)
      FCFlux = (/mixCFluxS(i1,j1),mixCFluxS(i1,j1p1),
     |    mixCFluxS(i1-1,j1),mixCFluxS(i1-1,j1p1) /)

      ting_pot_interp(i,j) = dot_product(mapping,FPot)
      ting_eng_interp(i,j) = max(dot_product(mapping,FEng),0.0)
      ting_flux_interp(i,j) = max(dot_product(mapping,FFlux),0.0)
      ting_bbeAvgEng_interp(i,j) = max(dot_product(mapping,FBEng),0.0)
      ting_bbeNumFlux_interp(i,j) =
     |   max(dot_product(mapping,FBFlux),0.0)
      ting_cuspAvgEng_interp(i,j) = max(dot_product(mapping,FCEng),0.0)
      ting_cuspNumFlux_interp(i,j) =
     |   max(dot_product(mapping,FCFlux),0.0)
      if (mytid==0.and.
     |   abs(ting_pot_interp(i,j)).gt.maxval(abs(FPot))*2.0) then
         write(6,*) "TGCM: OOB POTS: ", ting_pot_interp(i,j),i,j,i1,j1
         write(6,*) "TGCM: OOB LOCS: ", tempLon,tempLat,mixLon(j1,i1),
     |    mixLat(j1,i1)
         write(6,*) "TGCM: OOB VAL: ", FPOT
      endif
            !endif
         enddo
      enddo
      if (mytid == 0) then
      write(6,*) "TGCM: mix: min/max potS: ",
     |    minval(ting_pot_interp),maxval(ting_pot_interp)
      endif

      ! longitude pole
      ting_pot_interp(:,nlonp1) = ting_pot_interp(:,1)
      ting_eng_interp(:,nlonp1) = ting_eng_interp(:,1)
      ting_flux_interp(:,nlonp1) = ting_flux_interp(:,1)
      ting_bbeAvgEng_interp(:,nlonp1) = ting_bbeAvgEng_interp(:,1)
      ting_bbeNumFlux_interp(:,nlonp1) = ting_bbeNumFlux_interp(:,1)
      ting_cuspAvgEng_interp(:,nlonp1) = ting_cuspAvgEng_interp(:,1)
      ting_cuspNumFlux_interp(:,nlonp1) = ting_cuspNumFlux_interp(:,1)

      ! convert to correct energy
      ting_eng_interp = ting_eng_interp *0.5 !mean energy -> characteristic energy
      ting_bbeAvgEng_interp = ting_bbeAvgEng_interp*0.5
      ting_cuspAvgEng_interp = ting_cuspAvgEng_interp*0.5

      !call closeGroup()
      !call closeMixfile()
      write(6,*) "TGCM: Closing file, done reading"

      !write(6,*) "TGCM: Deleting ",mix2gcmLock
      !call CheckAndKill(mix2gcmLock)
!
! Process the imported data.
!

      endif ! mytid

      !write(*,*) " I AM BCAST: ",mytid

      call mpi_ibcast(ting_pot_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(1),ierr)
      call mpi_ibcast(ting_eng_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(2),ierr)
      call mpi_ibcast(ting_flux_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(3),ierr)
      call mpi_ibcast(ting_bbeAvgEng_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(4),ierr)
      call mpi_ibcast(ting_bbeNumFlux_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(5),ierr)
      call mpi_ibcast(ting_cuspAvgEng_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(6),ierr)
      call mpi_ibcast(ting_cuspNumFlux_interp,nlat*nlonp1,MPI_REAL8,0,
     |     TIEGCM_WORLD,req(7),ierr)
      call mpi_waitall(7,req,MPI_STATUSES_IGNORE,ierr)

      !write(*,*) " I AM DONE BCAST: ",mytid

      gpot = 0.
      geng = 0.
      gflx = 0.
      gpotm = 0.

      do jj=1,nlat
         do ii=1,nlon
            ! Note ii+2 skips over period points.
            gpot(ii+2,jj)=ting_pot_interp(jj,ii) !/1.2
            geng(ii+2,jj)=ting_eng_interp(jj,ii)
            gflx(ii+2,jj)=ting_flux_interp(jj,ii)
!            gBbeAvgEng(ii+2,jj)=ting_bbeAvgEng_interp(jj,ii)
!            gBbeNumFlux(ii+2,jj)=ting_bbeNumFlux_interp(jj,ii)
!            gCuspAvgEng(ii+2,jj)=ting_cuspAvgEng_interp(jj,ii)
!            gCuspNumFlux(ii+2,jj)=ting_cuspNumFlux_interp(jj,ii)
         enddo
         do ii=1,nlonp1
            gpotm(ii+2,jj)=ting_pot_interp(jj,ii)
         enddo
      enddo

      aurllbN = maxval(mixaurorabc1d(:,1))
      aurllbS = maxval(mixaurorabc1d(:,2))

      end subroutine ca_receive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_close
      use mpi_module,only: mytid
!
! This subroutine close cism ad hoc coupling interface for data exchange, and
! close Fortran I/O units
!
      if (mytid .eq. 0) then
        write(6,*)'finished closing Fortran Units for coupling'

        call closeMixFile
      endif

      end subroutine ca_close

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Function to check file existence and open NetCDF file
      subroutine openMixFile(ncFile, gStr)
          character(len=*), intent(in) :: ncFile
          character(len=*), optional :: gStr
          integer :: ierr, fExist

          ! Check if the file exists
          inquire(file=ncFile, exist=fExist)
          if (.not. fExist) then
              write(6, *) 'Unable to open file: ', trim(ncFile)
              stop
          endif

          ! Open the NetCDF file
          ierr = nf90_open(ncFile, NF90_NOWRITE, h5fId)
          write(*,*) "OPEN FILE: ",ierr

          ! Set initial group ID
          inId = h5fId

          ! If a group is specified, open it
          if (present(gStr)) then
              call openGroup(gStr)
          endif
      end subroutine openMixFile

      ! Function to open NetCDF group
      subroutine openGroup(gStr)
          character(len=*), intent(in) :: gStr
          integer :: ierr, gId

          ! Open the specified group
          !ierr = nf90_inq_ncid(ncid, trim(gStr), gId)
          write(*,*) "OPENING GROUP: ",trim(gStr),h5fId
          ierr = nf90_inq_grp_ncid(h5fId, trim(gStr), gId)
          write(*,*) "OPENED GROUP: ",trim(gStr),gId,h5fId,ierr
          if (ierr == nf90_NOERR) write(*,*) "NOERR"
          if (ierr == nf90_EBADID) write(*,*) "nf90_EBADID"
          if (ierr == nf90_ENOTNC4) write(*,*) "nf90_ENOTNC4"
          if (ierr == nf90_ESTRICTNC3) write(*,*) "nf90_ESTRICTNC3"
          if (ierr == nf90_EHDFERR) write(*,*) "nf90_EHDFERR"
          inId = gId
      end subroutine openGroup

      subroutine closeGroup()
          ! Reassign inId back to the root Id
          inId = h5fId
      end subroutine closeGroup

      ! Function to read 2D variable from NetCDF file
      subroutine read2DVariable(idStr, Q, gStrO)
          logical :: dsExist
          integer :: herr, dsTest
          integer :: typeClass, Nr, N
          integer :: typeSize
          !integer, allocatable, dimension(:) :: dims
          integer :: dims(2),dimsId(2)
          character(len=*), intent(in) :: idStr
          real*8, dimension(:,:), allocatable :: readArray
          real*8, dimension(:,:), allocatable :: Q, tmp
          character(len=*), intent(in), optional :: gStrO
          integer :: dsId,ierr
          character(len=mxlen_filename) :: tmpStr

          ! Check if this variable is a dataset or attribute
          ! Use nc_inq_varid to check the existence of the variable
          !write(*,*) "2D: ",trim(idStr),inId
          dsTest = nf90_inq_varid(inId, trim(idStr), dsId)
          write(*,*) "2DREAD: ",trim(idStr),inId,dsId
          if (dsTest == nf90_NOERR) write(*,*) "NOERR"
          if (dsTest == nf90_EBADID) write(*,*) "nf90_EBADID"
          if (dsTest == nf90_ENOTNC4) write(*,*) "nf90_ENOTNC4"
          if (dsTest == nf90_ESTRICTNC3) write(*,*) "nf90_ESTRICTNC3"
          if (dsTest == nf90_EHDFERR) write(*,*) "nf90_EHDFERR"
          dsExist = (dsTest == nf90_NOERR)

          if (dsExist) then
              ! Dataset exists, proceed with reading
              ierr = nf90_inquire_variable(inId, dsId, tmpStr,
     |           typeClass, Nr, dimsId, typeSize)
              ierr = nf90_inquire_dimension(inId,dimsId(1),
     |  tmpStr,dims(1))
              ierr = nf90_inquire_dimension(inId,dimsId(2),
     |  tmpStr,dims(2))

              !write(*,*) "TEST NAME: ",trim(tmpStr)
              !write(*,*) "TEST DIMS: ",dims
              !N = product(dims)

              ! Set relevant values in input chain and allocate memory
              if (allocated(readArray)) then
                  deallocate(readArray)
              endif
              !allocate(readArray(N))
              if (.not. allocated(readArray)) then
                allocate(readArray(dims(1),dims(2)))
              endif

              ! Read data
              ierr = nf90_get_var(inId, dsId, readArray)

              if (.not.allocated(Q)) then
                  write(*,*) "Allocate Q: "
                  allocate(Q(dims(2), dims(1)))
              endif
              !if (.not.allocated(tmp)) then
              !    allocate(tmp(dims(1), dims(2)))
              !endif

              !tmp = reshape(readArray, [dims(1), dims(2)])
              Q = transpose(readArray)
          else
              write(6, *) "Could not find this dataset: ",
     |    trim(idStr)
          endif
            !if (allocated(tmp)) then
            !    deallocate(tmp)
            !endif
            if (allocated(readArray)) then
                deallocate(readArray)
            endif
      end subroutine read2DVariable

      ! Function to close NetCDF file
      subroutine closeMixFile
          integer :: ierr

          ! Close the NetCDF file
          ierr = nf90_close(h5fId)
      end subroutine closeMixFile

      !subroutine readVariable(ncid,varname,varout)
      !character(len=*) :: varname
      !integer :: idv_var,istat,ncid
      !real*8,dimension(:,:) :: varout
      !
      !istat = nf90_inq_varid(ncid,varname,idv_var)
      !if (istat /= nf90_NOERR) call handle_ncerr(istat,
      !|  'ca_receive: Error getting '//trim(varname)//' var id')
      !varout = 0. ! init
      !istat = nf90_get_var_double(ncid,idv_var,varout)
      !if (istat /= nf90_NOERR) call handle_ncerr(istat,
      !|  'ca_receive: Error getting variable '//trim(varname))
      !
      !end subroutine readVariable
      !subroutine readVariable(varname,data)
      !character(len=*) :: varname
      !real, allocatable :: data(:,:)
      !real(4), allocatable :: readdata(:,:)
      !
      !integer(4) :: errorId
      !integer(4) :: varId
      !integer(4) :: selectId
      !integer(4) :: datatype
      !integer(4) :: nAttributes
      !integer(4) :: rank
      !integer(4) :: dimensions(32)
      !integer(4) :: indexStart(32)
      !integer(4) :: indexEnd(32)
      !character(len=1024) :: dataName
      !integer(4) :: i
      !
      !integer(4) sfn2index,sfselect,sfginfo
      !
      !integer(4) nPoints
      !
      !do i=1,32
      !   dimensions(i) = 0
      !   indexStart(i) = 0
      !   indexEnd(i) = 0
      !enddo
      !
      ! SDnametoindex
      !write(6,*) 'reading: ',varName
      !varId = sfn2index(fileId2, varName)
      !call errorCheck(varId,'sfn2index')
      ! sdselect
      !selectId = sfselect(fileId2, varId)
      !call errorCheck(selectId,'sfselect')
      ! sdgetinfo
      !errorId = sfginfo(selectId,dataName,rank, dimensions, datatype,nAttributes)
      !call errorCheck(errorId,'sfginfo')
      !
      !if (rank .ne. 2) then
      !   STOP "Only works with 2d variables!"
      !endif
      !do i=0,rank
      !   indexStart(i+1) = 0
      !   indexEnd(i+1) = dimensions(i+1)
      !enddo
      !
      ! Allocate memory
      !if (.not. allocated(data)) then
      !   allocate(data (dimensions(1), dimensions(2) ))
      !endif
      !if (.not. allocated(readdata)) then
      !   allocate(readdata (dimensions(1), dimensions(2) ))
      !endif
      !
      !errorId = sfrdata(selectId,indexStart,(/1,1/),indexEnd,readdata)
      !call errorCheck(errorId, "sfrdata")
      !
      !errorId = sfendacc(selectId)
      !call errorCheck(errorId, 'sfendacc')
      !
      !write(6,*) "On read: ", varName,minval(readdata),maxval(readdata)
      !data = real(readdata)
      !write(6,*) "On read2: ", varName, minval(data),maxval(data)
      !
      !end subroutine readVariable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !subroutine checkExistence(varName,exists)
      !character(len=*) :: varName
      !logical exists !output
      !integer(4) :: varId
      !
      !integer(4) sfn2index
      !
      ! SDnametoindex
      !varId = sfn2index(fileId2, varName)
      !if (varId.eq.-1) then
      !   exists=.false.
      !else
      !   exists=.true.
      !endif
      !end subroutine checkExistence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mix_search(p1,t1,p2,t2,j,i)
        ! for a given point (p2,t2) finds (j,i) of the grid (p1,t1)
        real*8,dimension(:,:),intent(in) :: p1,t1
        real*8,intent(in) :: p2,t2
        integer, intent(out) :: i,j
        integer :: Np

        Np = nlon

        ! FIXME: treat when outside the boundary
        j = minloc(abs(p1(:,2)-p2),1)
        i = minloc(abs(t1(2,:)-t2),1)

        ! correct for when we're close but just below this makes a big
        ! difference -- smoothes out the interpolation and removes
        ! spikes
        if ( (p1(j,2)-p2).gt.0 ) j=j-1
        if ( (t1(2,i)-t2).gt.0 ) i=i-1

        ! treat periodic
        if ( ((p2-p1(Np,2)).gt.0).and.(p2.lt.2*pi) ) j=Np
      end subroutine mix_search

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine tgcm_search(p1,t1,p2,t2,j,i)
        ! p is longitude
        ! t is latitude
        real*8,dimension(:,:),intent(in) :: p1,t1
        real*8,intent(in) :: p2,t2
        integer, intent(out) :: i,j
        integer :: Np,Nt

        Np = nmixlon
        Nt = nmixlat
        j = minloc(abs(p1(:,2)-p2),1)
        i = minloc(abs(t1(2,:)-t2),1)

        ! correction for when we're close but just below the point
        if ( (p1(j,2)-p2).gt. 0. ) j=j-1
        if ( (t1(2,i)-t2).gt. 0. ) i=i+1

        ! treat periodic
        if ( ((p2-p1(Np,1).gt.0).and.(p2.lt.2*pi) ) ) j=Np
        if (i .eq. 1) then
          write(*,*) "TOP?: ",t1(j,i),t2
          i = 2
        endif

      end subroutine tgcm_search

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine matinv4(A,B)
        !! Performs a direct calculation of the inverse of a 4×4 matrix.
        real*8, intent(in) :: A(4,4)   !! Matrix
        real*8, intent(out):: B(4,4)   !! Inverse matrix
        real*8             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv =
     |     1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+
     |       A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+
     |     A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
     |      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+
     |       A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+
     |     A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
     |      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+
     |       A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+
     |     A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
     |      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+
     |       A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+
     |     A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        ! Calculate the inverse of the matrix
        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+
     |                   A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+
     |                   A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+
     |                   A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+
     |                   A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+
     |                   A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+
     |                   A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+
     |                   A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+
     |                   A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+
     |                   A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+
     |                   A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+
     |                   A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+
     |                   A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+
     |                   A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+
     |                   A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+
     |                   A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+
     |                   A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+
     |                   A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+
     |                   A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+
     |                   A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+
     |                   A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+
     |                   A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+
     |                   A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+
     |                   A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+
     |                   A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+
     |                   A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+
     |                   A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+
     |                   A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+
     |                   A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+
     |                   A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+
     |                   A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+
     |                   A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+
     |                   A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
      end subroutine matinv4
!-----------------------------------------------------------------------

      subroutine CheckFileOrDie(fIn,errStr)
        character(len=*), intent(in) :: fIn,errStr

        logical :: fExist
        inquire(file=trim(fIn),exist=fExist)
        if (.not. fExist) then
            write(*,*) trim(errStr)
            write(*,*) "File: ", trim(fIn)
            write(*,*) ''
            stop
        endif

      end subroutine CheckFileOrDie

      !Get number of groups of form "Step#XXX" and start/end
      subroutine StepInfo(fStr, s0, sE, Nstp)
          use mpi_module,only : mytid,TIEGCM_WORLD
          use mpi
          character(len=*), intent(in) :: fStr
          integer, intent(out) :: s0, sE, Nstp
          integer ::  ierr, i, gExist
          logical :: isEnd, idFirst
          character(len=240) :: gStr
          ! This is approximately one carrington rotations at 1 minute outputs
          ! This should be on the upper limit of a magnetosphere run, but who knows!
          integer, parameter :: MAXSTEP0 = 40000

          isEnd = .false.
          idFirst = .false.
          i = 1

          ! Check if the file exists
          call CheckFileOrDie(fStr, "Unable to open file")
          !write(*,*) "Opening file",mytid,trim(fStr)
          ! Open the NetCDF file
          !ierr = nf90_open(trim(fStr), NF90_NOWRITE, ncid)
          !write(*,*) "Done open file",mytid

          do while (.not. isEnd)
              gExist = 0
              !write(*,*) "Loop step: ",mytid,i
              if (step1140Fix .and. i == 1140) then
                i = i+1
                cycle
              endif
              write(gStr, '(A,I0)') "Step#", i
              !write(6, '(A,I0)') "Step#", i
              ! Check if group exists in NetCDF file
              ierr = nf90_inq_grp_ncid(h5fId, trim(gStr), gExist)
              !write(*,*) "Group Found: ",ierr,NF90_NOERR,NF90_ENOGRP

              if (ierr .ne. NF90_NOERR) then
                  ! Step doesn't exist, either haven't found first or passed last
                  if (idFirst) then
                      ! Already found first, so we passed the last
                      !write(*,*) "FIRST GROUP: ",i,idFirst,isEnd
                      isEnd = .true.
                  else
                      ! Haven't yet found first, keep going
                      i = i + 1
                  endif
              else
                  ! Group exists, check if first
                  if (.not. idFirst) then
                      s0 = i
                      idFirst = .true.
                  endif
                  ! Increment no matter what
                  i = i + 1
              endif

              ! Give up after a while
              if (i >= MAXSTEP0) then
                  write(*, *) "Error, can't read file ", fStr
                  write(*, *) "Bailing out"
                  stop
              endif
          enddo

          write(*,*) "StepInfo Found: ",s0,i

          ! Close the NetCDF file
          !ierr = nf90_close(ncid)

          sE = i - 1
          Nstp = sE - s0 + 1

      end subroutine StepInfo

      !Same as StepTimes but for MJDs
      subroutine StepMJDs(ncFile, s0, sE, MJDs)
          character(len=*), intent(in) :: ncFile
          integer, intent(in) :: s0, sE
          real*8, intent(inout) :: MJDs(1:sE-s0+1)
          integer :: n, Nstp, ierr
          character(len=240) :: gStr
          integer :: gId
          logical :: aX
          real*8 :: M

          MJDs = 0.0

          Nstp = sE - s0 + 1

          ! Open the NetCDF file
          !call openMixFile(ncFile)

          do n = 1, Nstp
              if (step1140Fix .and. n == 1140) then
                MJDs(n) = MJDs(n-1)
                cycle
              endif
              write(gStr, '(A,I0)') "Step#", s0 + n - 1
              call openGroup(trim(gStr))

              ! Check if attribute "MJD" exists
              !ierr = nf90_inq_attid(inId, NF90_GLOBAL, "MJD", gId)
              aX = .true.
              if (aX) then
                  ! Read the "MJD" attribute
                  ierr = nf90_get_att(inId, NF90_GLOBAL, "MJD", MJDs(n))
              else
                  MJDs(n) = -TINY
              endif

              write(*,*) "MJD: ",MJDs(n)

              ! Close the current NetCDF group
              !call closeGroup()
          enddo

          ! Close the NetCDF file
          !call closeMixFile()
      end subroutine StepMJDs

      end module mage_oneway
#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine oneway_null
      end subroutine oneway_null
!-----------------------------------------------------------------------
#endif
