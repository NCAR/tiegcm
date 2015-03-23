!      
      module hme_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!      
! Read HME boundary perturbations T, U, V, W, RHO and interpolate
!   to model grid and time, and calculates Z.
!
! HME tidal components: daily files with hourly data on the tiegcm dres grid
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
      integer :: ncid=0    ! netcdf file id of currently open wacxg5 data file
      integer :: ntime=0   ! number of histories on open wacxg5 file
!
      integer,parameter :: mxtime=24 ! max number of times per file (hourly and daily files)
! Year, day of year, and hour from wacxg5 file:
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
     |  hme_v, 
     |  hme_w, 
     |  hme_rho 
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
     |  tlb0, ulb0, vlb0, wlb0, rlb0,  ! data at time index itime0
     |  tlb1, ulb1, vlb1, wlb1, rlb1   ! data at time index itime1

  100 continue
      if (ncid==0) call readfile(ncfiles(ifile),ifile,ncid)

      yfrac  = yfrac_tgcm(iyear,modeltime)
      if (itime0 < itime1) yfrac0 = yfrac_hme(itime0,modeltime(1))
      yfrac1 = yfrac_hme(itime1,modeltime(1))
!      write(6,*) 'in modeltime_tgcm',iyear,modeltime
!      write(6,*) 'in modeltime_wacx',wacxg5_date(itime1),
!     |  wacxg5_datesec(itime1)/60./60.
!      write(6,'(a28,3(x,f15.10))') 'in get_wacxg5_hourly',yfrac,
!     |  yfrac0,yfrac1 
     
      if (yfrac >= yfrac0 .and. yfrac < yfrac1) then
        if (.not.isread(itime0)) then
          call read_hme(ncid,ifile,itime0,tlb0,ulb0,vlb0,wlb0, rlb0)  ! read in LB values
          isread(itime0) = .true.
        endif
        if (.not.isread(itime1)) then
          call read_hme(ncid,ifile,itime1,tlb1,ulb1,vlb1,wlb1,rlb1)  ! read in LB values
          isread(itime1) = .true.
        endif   
     
! Interpolate LB values    
        call interpdata_lbc(modeltime,yfrac,yfrac0,yfrac1,
     |      tlb0,ulb0,vlb0,wlb0,rlb0,tlb1,ulb1,vlb1,wlb1,rlb1)
!	
        hme_ncfile_curr = ncfiles(ifile) ! for history file
!
	
      else ! next timestep on file
        tlb0 = tlb1 ; ulb0 = ulb1 ; vlb0 = vlb1 ; wlb0 = wlb1
        rlb0 = rlb1
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
      integer :: istart1(1),icount1(1)
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
      allocate(hme_modeltime(ntime,3),stat=istat)
      istart1(1) = 1
      icount1(1) = ntime
!
! year:
      istat = nf_inq_varid(ncid,"year",idv_year)
      istat = nf_get_vara_int(ncid,idv_year,istart1,icount1,
     |  hme_year)
!
! modeltime(doy,hour,min):
      istat = nf_inq_varid(ncid,"modeltime",idv_modeltime)
      istat = nf_get_vara_int(ncid,idv_modeltime,(/istart1,0/),
     |  (/icount1,3/),hme_modeltime)

!     write(6,"('readfile: file=',a,' date=',/,(6i10))") 
!    |  trim(ncfile),hme_year
!     write(6,"('readfile: modeltime=',/,(6i10))") hme_modeltime

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
      isecs = hme_modeltime(it,2)*3600+hme_modeltime(it,3)*60+
     |    hme_modeltime(it,4)
      yfrac_hme=float(hme_year(it))+
     |  (float(hme_modeltime(it,1))-1+float(isecs)/(3600.*24.))/maxday
!
      end function yfrac_hme
!--------------------------------------------------------------------------
      subroutine read_hme(ncid,ifile,itime,t,u,v,w,rho)
      implicit none
!
! Read 2d data t,u,v,w, rho from hme file.
! THME [deg K]; UHME,VHME,WHME [m/s]; RhoHME [%] 
! The file is already open under ncid.
!
! Args:
      integer,intent(in) :: ncid,ifile,itime
      real,dimension(nlon,nlat),intent(out) :: t,u,v,w,rho
!
! Local:
      integer :: istat,id,start3(3),count3(3)

      start3(1:2) = 1
      start3(3) = itime
      count3(1) = nlon
      count3(2) = nlat
      count3(3) = 1
!      
      istat = nf_inq_varid(ncid,"THME",id)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting THME var id')
      istat = nf_get_vara_double(ncid,id,start3,count3,t)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting variable THME')
      istat = nf_inq_varid(ncid,"UHME",id)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting UHME var id')
      istat = nf_get_vara_double(ncid,id,start3,count3,u)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting variable UHME')
      istat = nf_inq_varid(ncid,"VHME",id)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting VHME var id')
      istat = nf_get_vara_double(ncid,id,start3,count3,v)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting variable VHME')
      istat = nf_inq_varid(ncid,"WHME",id)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting WHME var id')
      istat = nf_get_vara_double(ncid,id,start3,count3,w)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting variable WHME')
      istat = nf_inq_varid(ncid,"RhoHME",id)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting RhoHME var id')
      istat = nf_get_vara_double(ncid,id,start3,count3,rho)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'read_hme: Error getting variable RhoHME')
      
!     write(6,"('Read data: ifile ',i3,': itime',i1,'=',i3,
!    |  ' min,max T=',2e12.4,' U=',2e12.4,' V=',2e12.4,' Z=',2e12.4)")
!    |  ifile,index,itime,minval(t),maxval(t),minval(u),maxval(u),
!    |  minval(v),maxval(v),minval(z),maxval(z)

      end subroutine read_hme
!-----------------------------------------------------------------------
      subroutine interpdata_lbc(modeltime,yfrac,yfrac0,yfrac1,
     |  t0,u0,v0,w0,rho0,t1,u1,v1,w1,rho1)
!
! Interpolate to yfrac between t0,u0,v0,w0,rho0 at yfrac0 and t1,u1,v1,w1,rho1 at yfrac1.
! Return module data hme_t,hme_u,hme_v,hme_w,hme_rho w/ unit conversion where necessary.
! Function time_interp is in util.F.
!
      use mpi_module,only: lat0,lat1,lon0,lon1
      implicit none
!
! Args:
      integer,intent(in) :: modeltime(4)
      real,intent(in) :: yfrac,yfrac0,yfrac1
      real,dimension(nlon,nlat),intent(in) :: 
     | t0,u0,v0,w0,rho0,t1,u1,v1,w1,rho1
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
            hme_w(i,j)= time_interp(w0(i-2,j),w1(i-2,j),yfrac0,
     |        yfrac1,yfrac)
            hme_rho(i,j)= time_interp(rho0(i-2,j),rho1(i-2,j),yfrac0,
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
          hme_w(1:2,j) = hme_w(nlon+1:nlon+2,j)
          hme_w(nlon+3:nlon+4,j) = hme_w(3:4,j)
          hme_rho(1:2,j) = hme_rho(nlon+1:nlon+2,j)
          hme_rho(nlon+3:nlon+4,j) = hme_rho(3:4,j)
        enddo ! j=1,nlat
!
! Convert from m/s to cm/s and m to cm:
      hme_rho = hme_rho*0.01 	! convert from [%] to [-]
      hme_u = hme_u*100.
      hme_v = hme_v*100.
!
      end subroutine interpdata_lbc
!-----------------------------------------------------------------------
      end module hme_module      
!-----------------------------------------------------------------------
      subroutine cal_hme_z(istep,itp,tbgrd,zprt,tprt,rhoprt)
      
      use mpi_module,only: lon0,lon1,lat0,lat1
      use params_module,only: nlonp4,nlat
      use cons_module,only: 
     |   grav,           ! 870. [cm/s2] accel due to gravity
     |   boltz,          ! 1.38E-16 [ergs/K] boltzman's constant 
     |   rmassinv,
     !   avo             ! 6.023e23 [#/mol] Avogadro constant 
      use fields_module,only: o2,o1,barm
      implicit none
!
! calculates the zpert 
!   zprt = - R/g0*tbgrd*(pprt/p0 + tprt/tbgrd)
! with R = kb/m [cm2/s2]  kb Boltzman constant; m mass
!
      integer,intent(in):: istep,  ! step
     |          itp               ! index of previous step
      real,intent(in) :: 
     |    tbgrd(nlonp4,nlat), ! t_lbc
     |    tprt(nlonp4,nlat),  ! tn pert HME [k]
     |    rhoprt(nlonp4,nlat) ! rho_pertr/rho_bgrd  HME [-]
      real,intent(out):: zprt(nlonp4,nlat)
      
      integer :: i,j,k
      real :: barm_lb(lon0:lon1,lat0:lat1),
     |     barm_tmp(2)
      
!     if(istep == 1) then
! 
! barm = mean molecular weight (k+1/2): [g/mol]
! barm_lb [g]
!
        do j=lat0,lat1
          do i=lon0,lon1
            do k=1,2
              barm_tmp(k) = 1./
     |  	(o2(k,i,j,itp)*rmassinv(1)+o1(k,i,j,itp)*rmassinv(2)+
     |  	(1.-o2(k,i,j,itp)-o1(k,i,j,itp))*rmassinv(3))
            enddo
            barm_lb(i,j) = 1.5*barm_tmp(1)-0.5*barm_tmp(2)
            barm_lb(i,j) = barm_lb(i,j)/avo
          enddo
        enddo
!     else  ! not first timestep
!
! barm1 = barm(k=0) (linear extrapolation)
!
!       do j=lat0,lat1
!         do i=lon0,lon1
!           barm_lb(i,j) = 1.5*barm(1,i,j,itp)-0.5*barm(2,i,j,itp)
!           barm_lb(i,j) = barm_lb(i,j)/avo
!         enddo
!       enddo
!     endif ! end first timestep
!    
! kb/m/g*Tn  [ergs/K /g / (cm/s2) *K] = [cm]
      do j=lat0,lat1
        do i=lon0,lon1
          zprt(i,j) = -boltz/grav/barm_lb(i,j)*tbgrd(i,j)*
     |         (rhoprt(i,j)+tprt(i,j)/tbgrd(i,j))

!         write(6,"('cal_hme_z: istep=',i3,' j=',i3,' i=',i3,
!    |      ' barm_lb=',e12.4,' tprt=',e12.4,' tbgrd=',e12.4,
!    |      ' rhoprt=',e12.4,' zprt=',e12.4)") istep,j,i,
!    |      barm_lb(i,j),tprt(i,j),tbgrd(i,j),rhoprt(i,j),zprt(i,j)

!          write(6,*) 'HME1',i,j,boltz/grav/barm_lb(i,j)*tprt(i,j)
!          write(6,*) 'HME2',i,j,
!     |  	boltz/grav/barm_lb(i,j)*tbgrd(i,j)*
!     |         (rhoprt(i,j))
!          write(6,*) 'HME3',i,j,
!     |  	boltz/grav/barm_lb(i,j)*tbgrd(i,j)*
!     |         (tprt(i,j)/tbgrd(i,j))
        enddo
      enddo
!     write(6,"('cal_hme_z: istep=',i3,' barm_lb min,max=',2e12.4,
!    |  ' zprt min,max=',2e12.4)") istep,minval(barm_lb),
!    |  maxval(barm_lb),minval(zprt),maxval(zprt)
        
      end subroutine cal_hme_z
