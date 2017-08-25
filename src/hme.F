!      
      module hme_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!      
! Read HME boundary perturbations T, U, V, W, Z and interpolate
!   to model grid and time, and calculates Z.
!
! HME tidal components: daily files with hourly data (0UT-23UT)on the tiegcm dres grid
!
      use params_module,only: nlon,nlat,nlonp4,nlonp2
      use mpi_module,only: lon0,lon1,lat0,lat1
      use nchist_module,only:nc_open,nc_close,handle_ncerr
      use addfld_module,only: addfld
      use input_module ,only: 
     |  nfiles =>hme_nfiles,     ! number of hme files
     |  ncfiles=>hme_ncfiles,    ! hme file names
     |  mxlen_filename,
     |  start_year               ! starting year (4-digit integer yyyy)
      use init_module,only: iyear,istep
!      
      implicit none
!      
      integer :: ncid=0    ! netcdf file id of currently open data file
      integer :: ntime=0   ! number of histories on open file
!
      integer,parameter :: mxtime=24 ! max number of times per file (hourly and daily files)
! Year, day of year, and hour from hme file:
      integer,allocatable,dimension(:),save :: 
     |  hme_year        ! year
      integer,allocatable,dimension(:,:),save :: 
     |  hme_modeltime   ! (doy,hr,min)    
!
#include <netcdf.inc>

!
      real,dimension(nlonp4,nlat),save ::
     |  hme_z, 
     |  hme_t, 
     |  hme_u, 
     |  hme_v
!
! Current hme_ncfile being used for nudging (for histories):
      character(len=mxlen_filename) :: hme_ncfile_curr
!
      contains 
!-----------------------------------------------------------------------
      subroutine get_hme(modeltime)
!
! Arg:
      integer,intent(in) :: modeltime(4)
!
! Local:
      integer,save :: ifile=1,itime0=1,itime1=2
      real,save :: yfrac,        ! year fraction at modeltime
     |             yfrac0,yfrac1 ! year fraction at data itime0,1
      logical,save :: isread(mxtime)=.false.
      integer :: istat
      real,dimension(nlon,nlat) ::       ! LB
     |  tlb0, ulb0, vlb0, zlb0,  ! data at time index itime0
     |  tlb1, ulb1, vlb1, zlb1   ! data at time index itime1

  100 continue
      if (ncid==0) call readfile(ncfiles(ifile),ifile,ncid)

      yfrac  = yfrac_tgcm(iyear,modeltime)
      if (itime0 < itime1) yfrac0 = yfrac_hme(itime0,modeltime(1))
      yfrac1 = yfrac_hme(itime1,modeltime(1))
!      write(6,*) 'in modeltime_tgcm',iyear,modeltime
      write(6,'(a28,3(x,f15.10))') 'in get_hme',yfrac,
     |  yfrac0,yfrac1 
     
      if (yfrac >= yfrac0 .and. yfrac < yfrac1) then
        if (.not.isread(itime0)) then
          call read_hme(ncid,ifile,itime0,tlb0,ulb0,vlb0,zlb0)  ! read in LB values
          isread(itime0) = .true.
        endif
        if (.not.isread(itime1)) then
          call read_hme(ncid,ifile,itime1,tlb1,ulb1,vlb1,zlb1)  ! read in LB values
          isread(itime1) = .true.
        endif   
     
! Interpolate LB values    
        call interpdata_lbc(modeltime,yfrac,yfrac0,yfrac1,
     |      tlb0,ulb0,vlb0,zlb0,tlb1,ulb1,vlb1,zlb1)
!	
        hme_ncfile_curr = ncfiles(ifile) ! for history file
!
	
      else ! next timestep on file
        tlb0 = tlb1 ; ulb0 = ulb1 ; vlb0 = vlb1 ; zlb0 = zlb1
!	
        if (itime0 > itime1) then
          itime0 = 1
          itime1 = 2
          isread(:) = .false.
        else
          itime0 = itime1
          itime1 = itime1+1
        endif
!
! End of current file. Retain t0 data at yfrac0 from this last time,
! reset itime1=1, and go back up to open new file. Note that now 
! itime0 > itime1 when interpolating between last time on previous 
! file and first time on new file, see conditionals above.
!
        if (itime1 > ntime) then
          ifile = ifile+1
          if (ifile > nfiles) stop 'out of hme files'
          istat = nf_close(ncid)
          ncid = 0 
          itime1 = 1 ! now itime0==ntime, so itime0 > itime1
          isread(itime1) = .false.
        endif
        goto 100
      endif        
!
      end subroutine get_hme
!-----------------------------------------------------------------------
      subroutine readfile(ncfile,ifile,ncid)
!
! Open hme file 
! Verify dimensions and coordinates, and read date and time
! but not LB values
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: ncfile
      integer,intent(in) :: ifile
      integer,intent(out) :: ncid
!
! Local:
      integer :: istat,ndims,idunlim,ngatts,nvars
      integer :: nlat_rd,nlon_rd
      integer :: id_lat,id_lon,id_time
      integer :: idv_year,idv_modeltime
      integer :: istart1(1),icount1(1),istart2(2),icount2(2)
      character(len=1024) :: dskfile
!
      dskfile = ' '
      call getfile(ncfile,dskfile)
      istat = nf_open(ncfile,NF_NOWRITE,ncid)
      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
      write(6,"('Opened hme file ',a)") trim(ncfile)
!
! Get and verify latitude dimension. 
! (read nlat_rd, compare with nlat, which is in params.h):
      istat = nf_inq_dimid(ncid,'lat',id_lat)
      istat = nf_inq_dimlen(ncid,id_lat,nlat_rd)
      if (nlat_rd /= nlat) then
        write(6,"(/,'>>> readfile: bad nlat_rd=',i3,
     |    ' -- should be nlat=',i3)") nlat_rd,nlat
        call shutdown('readfile')
      endif
!
! Get and verify longitude dimension. 
! (read nlon_rd, compare with nlon, which is in params.h):
      istat = nf_inq_dimid(ncid,'lon',id_lon)
      istat = nf_inq_dimlen(ncid,id_lon,nlon_rd)
      if (nlon_rd /= nlon) then
        write(6,"(/,'>>> readfile: bad nlon_rd=',i3,
     |    ' -- should be nlon=',i3)") nlon_rd,nlon
        call shutdown('readfile')
      endif  
!
! Get number of times (value of unlimited dimension) on the file:
      istat = nf_inq_unlimdim(ncid,id_time)     ! id of unlimited record var
      istat = nf_inq_dimlen(ncid,id_time,ntime) ! length of time var
      if (ntime > mxtime) then
        write(6,"('>>> Reading hme data file: ntime=',i4,' > ',
     |    'mxtime=',i4)") ntime,mxtime
        call shutdown('hme ntime')
      endif
!
! Allocate date and datesec vars:
      if (allocated(hme_year)) deallocate(hme_year)
      if (allocated(hme_modeltime)) deallocate(hme_modeltime)
      allocate(hme_year(ntime),stat=istat)
      allocate(hme_modeltime(3,ntime),stat=istat)
      istart1(1) = 1
      icount1(1) = ntime
!
! year:
      istat = nf_inq_varid(ncid,"year",idv_year)
      istat = nf_get_vara_int(ncid,idv_year,istart1,icount1,
     |  hme_year)
      if (istat /= NF_NOERR) then
        write(6,"('>>> Reading hme data file: year')") 
        call shutdown('hme year')
      endif
!
! modeltime(doy,hour,min):
      istart2(1) = 1
      istart2(2) = istart1(1)
      icount2(1) = 3
      icount2(2) = icount1(1)
      istat = nf_inq_varid(ncid,"mtime",idv_modeltime)
      istat = nf_get_vara_int(ncid,idv_modeltime,istart2,
     |  icount2,hme_modeltime)
      if (istat /= NF_NOERR) then
        write(6,"('>>> Reading hme data file: mtime')")
        call shutdown('hme mtime')
      endif

!     write(6,"('readfile: file=',a,' date=',/,(6i10))") 
!    |  trim(ncfile),hme_year
!      write(6,"('readfile: modeltime=',/,(6i10))") hme_modeltime(2,:)

      end subroutine readfile      
!-----------------------------------------------------------------------
      real function yfrac_tgcm(iyear,modeltime)
!
! Convert timegcm modeltime (day,hour,min) to fraction of a year. 
!
! Args:
      integer,intent(in) :: iyear,modeltime(4)
!
! Local:
      integer :: iday,isecs,iyear_loc
      real :: maxday
!
      iyear_loc = iyear
      iday      = modeltime(1)
!      
      if((mod(start_year,4).eq.0.and.mod(start_year,100).ne.0).or.
     |   (mod(start_year,400).eq.0)) then
         maxday = 366.                ! counting starts with 0
	 if(modeltime(1).eq.367) then ! should be day 1 of next year
	    iyear_loc = iyear-1
	 endif   
      else
         maxday = 365.  ! counting starts with 0
	 if(modeltime(1).eq.366) then ! should be day 1 of next year
	    iyear_loc = iyear-1
	 endif   
      endif
!
      isecs = modeltime(2)*3600+modeltime(3)*60+modeltime(4)
      yfrac_tgcm=float(iyear_loc)+
     |  (float(iday)-1+float(isecs)/(3600.*24.))/maxday
!
      end function yfrac_tgcm
!-----------------------------------------------------------------------
      real function yfrac_hme(it,itgcm_day)
!
! Convert hme date and time to fraction of a year. 
!
! Args:
      integer,intent(in) :: it,itgcm_day
!
! Local:
      integer :: iday,isecs
      real :: maxday
!      
      if((mod(start_year,4).eq.0.and.mod(start_year,100).ne.0).or.
     |   (mod(start_year,400).eq.0).and.(itgcm_day.ne.367)) then
         maxday = 366.                ! counting starts with 0	 
      else
         maxday = 365.  ! counting starts with 0
      endif
!
      isecs = hme_modeltime(2,it)*3600+hme_modeltime(3,it)*60
      yfrac_hme=float(hme_year(it))+
     |  (float(hme_modeltime(1,it))-1+float(isecs)/(3600.*24.))/maxday
!       write(6,*) 'yfrac_hme ',hme_modeltime(1,it),
!     |   hme_modeltime(2,it),hme_modeltime(3,it)
!
      end function yfrac_hme
!--------------------------------------------------------------------------
      subroutine read_hme(ncid,ifile,itime,t,u,v,z)
      implicit none
!
! Read 2d data t,u,v,w, rho from hme file.
! THME [deg K]; UHME,VHME,WHME [m/s]; ZHME [m]
! The file is already open under ncid.
!
! Args:
      integer,intent(in) :: ncid,ifile,itime
      real,dimension(nlon,nlat),intent(out) :: t,u,v,z
!
! Local:
      integer :: istat,id,start3(3),count3(3)

      start3(1:2) = 1
      start3(3) = itime
      count3(1) = nlon
      count3(2) = nlat
      count3(3) = 1
!      
      istat = nf_inq_varid(ncid,"TN",id)
      if (istat /= NF_NOERR) then 
         call handle_ncerr(istat,
     |  'read_hme: Error getting TN var id')
         call shutdown('read_hme TN')
      endif	 
      istat = nf_get_vara_double(ncid,id,start3,count3,t)
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |  'read_hme: Error getting variable TN')
         call shutdown('read_hme TN')
      endif	 
      istat = nf_inq_varid(ncid,"UN",id)
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |  'read_hme: Error getting UN var id')
         call shutdown('read_hme UN')
      endif	 
      istat = nf_get_vara_double(ncid,id,start3,count3,u)
      if (istat /= NF_NOERR) then 
        call handle_ncerr(istat,
     |  'read_hme: Error getting variable UN')
         call shutdown('read_hme UN')
      endif	 
      istat = nf_inq_varid(ncid,"VN",id)
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |  'read_hme: Error getting VN var id')
         call shutdown('read_hme VN')
      endif	 
      istat = nf_get_vara_double(ncid,id,start3,count3,v)
      if (istat /= NF_NOERR) then 
        call handle_ncerr(istat,
     |  'read_hme: Error getting variable VN')
         call shutdown('read_hme VN')
      endif	 
      istat = nf_inq_varid(ncid,"Z",id)
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |  'read_hme: Error getting Z var id')
         call shutdown('read_hme Z')
      endif	 
      istat = nf_get_vara_double(ncid,id,start3,count3,z)
      if (istat /= NF_NOERR) then
        call handle_ncerr(istat,
     |  'read_hme: Error getting variable Z')
         call shutdown('read_hme Z')
      endif	

      end subroutine read_hme
!-----------------------------------------------------------------------
      subroutine interpdata_lbc(modeltime,yfrac,yfrac0,yfrac1,
     |  t0,u0,v0,z0,t1,u1,v1,z1)
!
! Interpolate to yfrac between t0,u0,v0,z0 at yfrac0 and t1,u1,v1,z1 at yfrac1.
! Return module data hme_t,hme_u,hme_v,hme_z w/ unit conversion where necessary.
! Function time_interp is in util.F.
!
      use mpi_module,only: lat0,lat1,lon0,lon1
      implicit none
!
! Args:
      integer,intent(in) :: modeltime(4)
      real,intent(in) :: yfrac,yfrac0,yfrac1
      real,dimension(nlon,nlat),intent(in) :: 
     | t0,u0,v0,z0,t1,u1,v1,z1
!
! Local:
      integer :: i,j,k
! External:

      real,external :: time_interp ! util.F

        do j=1,nlat
          do i=3,nlon+2
            hme_t(i,j)= time_interp(t0(i-2,j),t1(i-2,j),yfrac0,
     |        yfrac1,yfrac)
            hme_u(i,j)= time_interp(u0(i-2,j),u1(i-2,j),yfrac0,
     |        yfrac1,yfrac)
            hme_v(i,j)= time_interp(v0(i-2,j),v1(i-2,j),yfrac0,
     |        yfrac1,yfrac)
            hme_z(i,j)= time_interp(z0(i-2,j),z1(i-2,j),yfrac0,
     |        yfrac1,yfrac)
          enddo ! i=1,nlon
!
! Periodic points:
          hme_t(1:2,j) = hme_t(nlon+1:nlon+2,j)
          hme_t(nlon+3:nlon+4,j) = hme_t(3:4,j)
          hme_u(1:2,j) = hme_u(nlon+1:nlon+2,j)
          hme_u(nlon+3:nlon+4,j) = hme_u(3:4,j)
          hme_v(1:2,j) = hme_v(nlon+1:nlon+2,j)
          hme_v(nlon+3:nlon+4,j) = hme_v(3:4,j)
          hme_z(1:2,j) = hme_z(nlon+1:nlon+2,j)
          hme_z(nlon+3:nlon+4,j) = hme_z(3:4,j)
        enddo ! j=1,nlat
!
! Convert from m/s to cm/s and m to cm:
      hme_u = hme_u*100. ! convert from m/s to cm/s
      hme_v = hme_v*100. ! convert from m/s to cm/s
      hme_z = hme_z*100. ! convert from m to cm
!
      end subroutine interpdata_lbc
!-----------------------------------------------------------------------
      end module hme_module 
