!      
      module bgrd_data_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!      
! Read background based on SABER,WINDII/HRDI at lower boundary Z, TN, UN, VN, 
!  and interpolate time.
!
! file name: saber_lbcs_tngph_windiihrdi_unvn_notides.nc
! only for dres (lon from -180,..177.5 by 2.5 deg; 
!                lat -88.75 to 88.75 by 2.5 deg)
!  double Z(time, nmonth, lat, lon) ;
!                Z:long_name = "Geopotential Height Perturbation" ;
!                Z:units = "meters" ;
!        double TN(time, nmonth, lat, lon) ;
!                TN:long_name = "Neutral Temperature Perturbation" ;
!                TN:units = "deg K" ;
!        double UN(time, nmonth, lat, lon) ;
!                UN:long_name = "Neutral Zonal Wind Perturbation" ;
!                UN:units = "m/s" ;
!        double VN(time, nmonth, lat, lon) ;
!                VN:long_name = "Neutral Meridional Wind Perturbation" ;
!                VN:units = "m/s" ;
!
      use params_module,only: nlon,nlat,nlonp4,nlonp2
      use mpi_module,only: lon0,lon1,lat0,lat1
      use nchist_module,only:nc_open,nc_close,handle_ncerr
      use addfld_module,only: addfld
!
      implicit none
!
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif 
! 
! Bachground based on data at lower boundary at model grid, 
! output by this module, to be used by dt, duv, addiag, etc.
! Allocated by sub alloc_bgrddata (called by allocdata):
!
      real,allocatable,dimension(:,:) :: ! (lon0:lon1,lat0:lat1)
     |  bgrddata_z,
     |  bgrddata_t,
     |  bgrddata_u,
     |  bgrddata_v
!
! Private module data, read by sub rdbgrddata:
! Will be allocated (lon0:lon1,lat0:lat1,nmonth,nhour)
      integer,parameter,private :: nmonth= 12, nhour = 24
      real,allocatable,dimension(:,:,:,:),private ::
     |  z_bgrd,
     |  t_bgrd,
     |  u_bgrd,
     |  v_bgrd
!
      contains
!-----------------------------------------------------------------------
      subroutine get_bgrddata(istep,iyear,iday,secs) 
!
! Module driver to read nc files, and do time interpolations.
!
! Files provided by user via namelist read:
      use input_module,only:  bgrddata_ncfile
!
! Integer flags set according to user-requested files:
      use init_module,only:  ibgrddata   ! 0/1 flag to get background data
!
! Driver for obtaining GSWM data, called once per timestep from advance.
!
! Args:
      integer,intent(in) :: istep,iyear,iday
      real,intent(in) :: secs
!
! Local:
      integer :: iprint
!
      iprint = 0
      if (istep==1) iprint = 1
!
! Get gswm migrating diurnal tide data:
      if (ibgrddata > 0) then
        if (istep==1) call rd_bgrddata(bgrddata_ncfile)
        call mk_bgrddata(iyear,iday,int(secs),iprint)
      endif
!
      end subroutine get_bgrddata 
! 
!-----------------------------------------------------------------------
      subroutine rd_bgrddata(ncfile)
      use input_module,only: mxlen_filename
      implicit none
!
! Args:
      character(len=*),intent(in) :: ncfile
!
! Local:
      integer :: ncid,istat
      character(len=mxlen_filename) :: dskfile
      integer :: nlon_rd, nlat_rd, nmonth_rd, nhour_rd
      integer :: id_nmonth, id_nhour, id_nlon, id_nlat
      integer :: idv_z, idv_t, idv_u, idv_v
      character(len=240) :: char240
      real,dimension(nlonp4,nlat,nmonth,nhour) :: t,u,v,z
      real,dimension(nlon,nlat,nmonth,nhour) :: tmp
!
      dskfile = ' '
      call getfile(ncfile,dskfile)
      write(6,"(/,72('-'))")
      write(6,"('Reading BGRDDATA data file ',a)") trim(ncfile)
!
! Open netcdf file:
      call nc_open(ncid,dskfile,'OLD','READ')
      if (ncid==0) then
        write(6,"(/,'>>> rd:_bgrddata error opening netcdf backgrund',
     |    ' data file ',a,' dskfile ',a)") trim(ncfile),trim(dskfile)
        call shutdown('rd_bgrddata')
      endif
!
! Check nmonth dimension:
      istat = nf_inq_dimid(ncid,'nmonth',id_nmonth)
      istat = nf_inq_dimlen(ncid,id_nmonth,nmonth_rd)
      if (istat /= NF_NOERR) then
        write(char240,"('rd_bgrddata: Error getting nmonth dimension',
     |    ' from file ',a)") trim(ncfile)
        call handle_ncerr(istat,char240)
      endif
      if (nmonth_rd /= nmonth) then
        write(6,"(/,'>>> rd_bgrddata: nmonth_rd=',i4,' not equal to',
     |	  ' nmonth=',i4)") nmonth_rd,nmonth
        write(6,"('background data file: ',a)") trim(ncfile)
        call shutdown('rd_bgrddata')
      endif
!
! Get nhour (time) dimension:
      istat = nf_inq_dimid(ncid,'time',id_nhour)
      istat = nf_inq_dimlen(ncid,id_nhour,nhour_rd)
      if (istat /= NF_NOERR) then
        write(char240,"('rd_bgrddata: Error getting time dimension',
     |    ' from file ',a)") trim(ncfile)
        call handle_ncerr(istat,char240)
      endif
      if (nhour_rd /= nhour) then
        write(6,"(/,'>>> rd_bgrddata: nhour_rd=',i4,' not equal to',
     |	  ' nhour=', i4)") nhour_rd,nhour
        write(6,"('background data file: ',a)") trim(ncfile)
        call shutdown('rd_bgrddata')
      endif
!
! Get nlon dimension:
      istat = nf_inq_dimid(ncid,'lon',id_nlon)
      istat = nf_inq_dimlen(ncid,id_nlon,nlon_rd)
      if (istat /= NF_NOERR) then
        write(char240,"('rd_bgrddata: Error getting nlon dimension',
     |    ' from file ',a)") trim(ncfile)
        call handle_ncerr(istat,char240)
      endif
      if (nlon_rd /= nlon) then
        write(6,"(/,'>>> rd_bgrddata: nlon_rd=',i4,' not equal to',
     |	  ' nlon=',i4)") nlon_rd,nlon
        write(6,"('background data file: ',a)") trim(ncfile)
        call shutdown('rd_bgrddata')
      endif
!
! Get nlat dimension:
      istat = nf_inq_dimid(ncid,'lat',id_nlat)
      istat = nf_inq_dimlen(ncid,id_nlat,nlat_rd)
      if (istat /= NF_NOERR) then
        write(char240,"('rd_bgrddata: Error getting nlat dimension',
     |    ' from file ',a)") trim(ncfile)
        call handle_ncerr(istat,char240)
      endif
      if (nlat_rd /= nlat) then
        write(6,"(/,'>>> rd_bgrddata: nlat_rd=',i4,' not equal to',
     |	' nlat=',i4)") nlat_rd,nlat
        write(6,"('background data file: ',a)") trim(ncfile)
        call shutdown('rd_bgrddata')
      endif
!
! Get Z geopotential height perturbation [m]:
!
      istat = nf_inq_varid(ncid,'Z',idv_z)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting Z var id')
      z = 0. ! init
      istat = nf_get_var_double(ncid,idv_z,tmp)
      z(3:nlon+2,:,:,:) = tmp(:,:,:,:)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting variable Z')
      z = z*100. 	! convert from m to cm
!     write(6,"('rdgswm ',a,': z min,max=',2e12.4)") 
!    |  type,minval(z),maxval(z)
!
! Get TN perturbation [deg K]:
!
      istat = nf_inq_varid(ncid,'TN',idv_t)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting TN var id')
      t = 0.
      istat = nf_get_var_double(ncid,idv_t,tmp)
      t(3:nlon+2,:,:,:) = tmp(:,:,:,:)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting variable TN')
!     write(6,"('rdgswm ',a,': tn min,max=',2e12.4)") 
!    |  type,minval(t),maxval(t)
!
! Get UN perturbation [m/s]:
!
      istat = nf_inq_varid(ncid,'UN',idv_u)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting UN var id')
      u = 0.
      istat = nf_get_var_double(ncid,idv_u,tmp)
      u(3:nlon+2,:,:,:) = tmp(:,:,:,:)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting variable UN')
      u = u*100. 	! convert to cm/s
!     write(6,"('rdgswm ',a,': un min,max=',2e12.4)") 
!    |  type,minval(u),maxval(u)
!
! Get VN perturbation [m/s]:
!
      istat = nf_inq_varid(ncid,'VN',idv_v)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting VN var id')
      v = 0.
      istat = nf_get_var_double(ncid,idv_v,tmp)
      v(3:nlon+2,:,:,:) = tmp(:,:,:,:)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |  'rd_bgrddata: Error getting variable VN')
      v = v*100. 	! convert to cm/s
!     write(6,"('rdgswm ',a,': vn min,max=',2e12.4)") 
!    |  type,minval(v),maxval(v)
!
! Transfer to private module data (whole-array assignments):
      z_bgrd(lon0:lon1,lat0:lat1,:,:) = z(lon0:lon1,lat0:lat1,:,:)
      t_bgrd(lon0:lon1,lat0:lat1,:,:) = t(lon0:lon1,lat0:lat1,:,:)
      u_bgrd(lon0:lon1,lat0:lat1,:,:) = u(lon0:lon1,lat0:lat1,:,:)
      v_bgrd(lon0:lon1,lat0:lat1,:,:) = v(lon0:lon1,lat0:lat1,:,:)
!
! Close the file:
      call nc_close(ncid)
      write(6,"('Completed read from Background data file ',a)") 
     |   trim(ncfile)
      write(6,"(/,72('-'))")
      end subroutine rd_bgrddata
!-----------------------------------------------------------------------
      subroutine mk_bgrddata(iyear,iday,isecs,iprint)
!
! Use data read from background data file to return  U,V, Z
! and T at current model date and time. It is assumed that the background file
! provides data for a whole year with one day of data for each month
! with hourly values (0 UT to 23 UT).
! The data is linearly interpolated first to the model UT for the
! corresponding months (current and next month) and then linearly
! interpolated to the modelday
!
      use hist_module,only: modeltime
#ifdef MPI
      use mpi_module,only: mp_periodic_f2d
#endif
      implicit none
!
! Args:
      integer,intent(in) :: iyear,iday,isecs,iprint
!
! Local:
      real :: difsec,difday,sec,secint,dayint
      integer :: i,j,ndmon_nl(13),ndmon_lp(13),nonleap,
     |   mon_cur,mon_nxt,ihr_cur,ihr_nxt,nointp=0,
     |   lonbeg,lonend,i0,i1,j0,j1
      real,dimension(lon0:lon1,lat0:lat1) :: t,u,v,z,
     |  z_curmo, t_curmo, u_curmo, v_curmo,
     |  z_nxtmo, t_nxtmo, u_nxtmo, v_nxtmo,
     |  z_curmo_nxthr, t_curmo_nxthr, u_curmo_nxthr, v_curmo_nxthr,
     |  z_nxtmo_nxthr, t_nxtmo_nxthr, u_nxtmo_nxthr, v_nxtmo_nxthr
      real :: fbgrd(lon0:lon1,lat0:lat1,4) ! for mp calls
!   
! Data: GSWM data given at the middle of each month 
!       -> day of year for this data assuming non-leap year
!         J  F  M  A   M   J   J   A   S   O   N   D   J
      data ndmon_nl                                ! non_leap
     |  /15,46,74,105,135,166,196,227,258,288,319,349,380/ 
! 
! External:
      real,external :: finterp
!
! For addfld calls:
      i0 = lon0 ; i1 = lon1
      j0 = lat0 ; j1 = lat1
! 
! Get month of model run
      do i = 1,13
        if(iday.le.ndmon_nl(i)) goto 10
      enddo
 10   mon_nxt = i		! next month
      if(mon_nxt == 13 ) mon_nxt = 1
      mon_cur = i-1		! current month
      if(mon_cur == 0 ) mon_cur = 12
! 
! Get hour of model run (model hours from 0 to 23 )
      ihr_cur = modeltime(2)	! current hour
      ihr_nxt = modeltime(2)+1  ! next hour
      if(ihr_nxt == 24 ) ihr_nxt = 0
!
! Subdomains at current month, current hour:
      z_curmo = bgrddat(mon_cur,ihr_cur+1,'z')
      t_curmo = bgrddat(mon_cur,ihr_cur+1,'t')
      u_curmo = bgrddat(mon_cur,ihr_cur+1,'u')
      v_curmo = bgrddat(mon_cur,ihr_cur+1,'v')
!
! Subdomains at next month, current hour:
      z_nxtmo = bgrddat(mon_nxt,ihr_cur+1,'z')
      t_nxtmo = bgrddat(mon_nxt,ihr_cur+1,'t')
      u_nxtmo = bgrddat(mon_nxt,ihr_cur+1,'u')
      v_nxtmo = bgrddat(mon_nxt,ihr_cur+1,'v')
!
! Interpolate to month:
      if (isecs /= 0) then
        difsec = 60.*60. ! difference in sec between ihr_cur and ihr_nxt
        secint = float(modeltime(3)*60 + modeltime(4))  ! interpolation time [sec]
!
! Current month, next hour:
        z = z_curmo ; t = t_curmo ; u = u_curmo ; v = v_curmo
        z_curmo_nxthr = bgrddat(mon_cur,ihr_nxt+1,'z')
        t_curmo_nxthr = bgrddat(mon_cur,ihr_nxt+1,'t')
        u_curmo_nxthr = bgrddat(mon_cur,ihr_nxt+1,'u')
        v_curmo_nxthr = bgrddat(mon_cur,ihr_nxt+1,'v')

        call timeinterp(z, z_curmo_nxthr, difsec, secint, z_curmo, 1)
        call timeinterp(t, t_curmo_nxthr, difsec, secint, t_curmo, 1)
        call timeinterp(u, u_curmo_nxthr, difsec, secint, u_curmo, 1)
        call timeinterp(v, v_curmo_nxthr, difsec, secint, v_curmo, 1)
!
! Interpolate to next month:
        z = z_nxtmo ; t = t_nxtmo ; u = u_nxtmo ; v = v_nxtmo
!
! Next month, next hour:
        z_nxtmo_nxthr = bgrddat(mon_nxt,ihr_nxt+1,'z')
        t_nxtmo_nxthr = bgrddat(mon_nxt,ihr_nxt+1,'t')
        u_nxtmo_nxthr = bgrddat(mon_nxt,ihr_nxt+1,'u')
        v_nxtmo_nxthr = bgrddat(mon_nxt,ihr_nxt+1,'v')

        call timeinterp(z, z_nxtmo_nxthr, difsec, secint, z_nxtmo,1)
        call timeinterp(t, t_nxtmo_nxthr, difsec, secint, t_nxtmo,1)
        call timeinterp(u, u_nxtmo_nxthr, difsec, secint, u_nxtmo,1)
        call timeinterp(v, v_nxtmo_nxthr, difsec, secint, v_nxtmo,1)

      endif ! interpolate to month or not
!
! Check if interpolation to ut is necessary
      nointp = 0
      if(iday.eq.ndmon_nl(mon_cur)) then ! same day as cur. month
         nointp = 1
         goto 20
      endif
      if(iday.eq.ndmon_nl(mon_nxt)) then ! same day as next month
         nointp = 2
         goto 20
      endif
!
! If ut interpolation is necessary, calculate time differences
      if(mon_cur /= 12) then                           ! not December
        difday = ndmon_nl(mon_nxt)-ndmon_nl(mon_cur)   ! difference in days
        dayint = iday - ndmon_nl(mon_cur) ! difference to interpolation day
      else                                             ! December wrap around
        difday = ndmon_nl(mon_cur+1)-ndmon_nl(mon_cur) ! difference in days
        if(iday.lt.ndmon_nl(mon_nxt)) then ! difference to interpolation day
           dayint = 365. - ndmon_nl(mon_cur)+  iday
        else
          dayint = iday - ndmon_nl(mon_cur)
        endif
      endif
 20   continue          ! if no interpolation is necessary (nointp /= 0)
!
! Interpolate to ut if necessary to local z,t,u,v:
      select case (nointp)
        case (0)                  ! interpolate
          call timeinterp(z_curmo,z_nxtmo,difday,dayint,z,1)
          call timeinterp(t_curmo,t_nxtmo,difday,dayint,t,1)
          call timeinterp(u_curmo,u_nxtmo,difday,dayint,u,1)
          call timeinterp(v_curmo,v_nxtmo,difday,dayint,v,1)
        case (1)                  ! no interp (same day as current month)
          z = z_curmo
          t = t_curmo
          u = u_curmo
          v = v_curmo
        case (2)                  ! no interp (same day as next month)
          z = z_nxtmo
          t = t_nxtmo
          u = u_nxtmo
          v = v_nxtmo
        case default
          write(6,"(/,'>>> mk_bgrddata: unknown nointp=',i4)") nointp
          call shutdown('mk_bgrddata')
      end select
!
! Transfer to module data according to type:
       bgrddata_z(lon0:lon1,lat0:lat1) = z
       bgrddata_t(lon0:lon1,lat0:lat1) = t
       bgrddata_u(lon0:lon1,lat0:lat1) = u
       bgrddata_v(lon0:lon1,lat0:lat1) = v
!
! Do mpi periodic points exchange for gswm with f2d(:)
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
#ifdef MPI
       fbgrd(:,:,1) = bgrddata_z(lon0:lon1,lat0:lat1)
       fbgrd(:,:,2) = bgrddata_t(lon0:lon1,lat0:lat1)
       fbgrd(:,:,3) = bgrddata_u(lon0:lon1,lat0:lat1)
       fbgrd(:,:,4) = bgrddata_v(lon0:lon1,lat0:lat1)
       call mp_periodic_f2d(fbgrd,lon0,lon1,lat0,lat1,4)
       bgrddata_z(lon0:lon1,lat0:lat1) = fbgrd(:,:,1)
       bgrddata_t(lon0:lon1,lat0:lat1) = fbgrd(:,:,2)
       bgrddata_u(lon0:lon1,lat0:lat1) = fbgrd(:,:,3)
       bgrddata_v(lon0:lon1,lat0:lat1) = fbgrd(:,:,4)
#else
       call set_periodic_f2d(bgrddata_z)
       call set_periodic_f2d(bgrddata_t)
       call set_periodic_f2d(bgrddata_u)
       call set_periodic_f2d(bgrddata_v)
#endif

       call addfld('bgrddata_z','bgrd Z','cm',
     |	 bgrddata_z(i0:i1,j0:j1),'lon',i0,i1,'lat',j0,j1,0)
       call addfld('bgrddata_t','bgrd TN','K',
     |	 bgrddata_t(i0:i1,j0:j1),'lon',i0,i1,'lat',j0,j1,0)
       call addfld('bgrddata_u','bgrd UN','cm/s',
     |	 bgrddata_u(i0:i1,j0:j1),'lon',i0,i1,'lat',j0,j1,0)
       call addfld('bgrddata_v','bgrd VN','cm/s',
     |	 bgrddata_v(i0:i1,j0:j1),'lon',i0,i1,'lat',j0,j1,0)

      end subroutine mk_bgrddata
!
!-----------------------------------------------------------------------
      function bgrddat(mon,ihr,fname)
!
! Return subdomain array from read arrays according to requested field, 
!   month, and hour.
! The 4 global fields t_bgrd, etc, are private module data, read by 
!   sub rd_bgrddata at beginning of model run.
!
! Args:
      character(len=*),intent(in) :: fname
      integer,intent(in) :: mon,ihr
!
! Function output dimension:
      real :: bgrddat(lon0:lon1,lat0:lat1)
!
! Field type must be t,u,v, or z:
      if (fname/='t'.and.fname/='u'.and.fname/='v'.and.fname/='z') then
        write(6,"(/,'>>> bgrddat: unknown fname=',a)") fname
        write(6,"('Must be t, u, v, or z')")
        call shutdown('bgrddat')
      endif
      if (fname=='t') bgrddat(:,:)= t_bgrd(lon0:lon1,lat0:lat1,mon,ihr)
      if (fname=='u') bgrddat(:,:)= u_bgrd(lon0:lon1,lat0:lat1,mon,ihr)
      if (fname=='v') bgrddat(:,:)= v_bgrd(lon0:lon1,lat0:lat1,mon,ihr)
      if (fname=='z') bgrddat(:,:)= z_bgrd(lon0:lon1,lat0:lat1,mon,ihr)
      
      end function bgrddat
!-----------------------------------------------------------------------
      subroutine timeinterp(d1,d2,difd1d2,difint,fout,
     |  iprnt)
!
! Interpolate from fields d1,d2  linearly to time difint
!   (d1 must be at 0 unit time), returning fout
! On input:
!   d1,d2   = the field (d1 at 0 unit, d2 at difd1d2 units)
!   difd1d2 = time difference [unit] between d1 & d2
!   difint  = time of interpolation [unit] (counted from d1=0 time)
! In output:
!   fout is defined at difint
!
! Args:
      integer,intent(in) :: iprnt
      real,intent(in)  :: difd1d2,difint,d1(lon0:lon1,lat0:lat1),
     |   d2(lon0:lon1,lat0:lat1)
      real,intent(out) :: fout(lon0:lon1,lat0:lat1)
!
! Local:
      integer :: i,j
      real :: frac0,difd1d2_inv
!
      fout = 0. ! initialize output      
!
! Interpolation:
! from data_cur/d1 to data_nxt/d2 the difference is in time difd1d2
!      
      difd1d2_inv = 1./difd1d2
      frac0 = difint*difd1d2_inv	! x(int)/[x(d2)-x(d1)]
 
! linear interpolation: special case x(d1) = 0    
!  fout = (d2-d1)*frac0 + d1 - (d2-d1)*difd1d2_inv*x(d1)
! (d2-d1)*difd1d2_inv*x(d1) = 0 since x(d1) = 0
      
      do i = lon0,lon1
        do j = lat0,lat1
          fout(i,j) = (d2(i,j)-d1(i,j))*frac0 + d1(i,j)
	enddo
      enddo
      
      end subroutine timeinterp
!-----------------------------------------------------------------------
      subroutine set_periodic_f2d(f)
!
! Set periodic points for all f2d fields (serial or non-mpi only):
!
      real,intent(inout) :: f(nlonp4,nlat)

! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
      f(1:2,:)              = f(nlonp4-3:nlonp4-2,:)
      f(nlonp4-1:nlonp4,:)  = f(3:4,:)
      end subroutine set_periodic_f2d
!-----------------------------------------------------------------------
      subroutine alloc_bgrddata(lon0,lon1,lat0,lat1)
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Module data to be use-associated to lbc.F:
      allocate(bgrddata_z (lon0:lon1,lat0:lat1),stat=istat)
      allocate(bgrddata_t (lon0:lon1,lat0:lat1),stat=istat)
      allocate(bgrddata_u (lon0:lon1,lat0:lat1),stat=istat)
      allocate(bgrddata_v (lon0:lon1,lat0:lat1),stat=istat)

      write(6,"('Allocated background t,u,v,z (lon0:lon1,lat0:lat1)')")
!
! Private module data to be read from nc file:
      allocate(z_bgrd (lon0:lon1,lat0:lat1,nmonth,nhour),stat=istat)
      allocate(t_bgrd (lon0:lon1,lat0:lat1,nmonth,nhour),stat=istat)
      allocate(u_bgrd (lon0:lon1,lat0:lat1,nmonth,nhour),stat=istat)
      allocate(v_bgrd (lon0:lon1,lat0:lat1,nmonth,nhour),stat=istat)

      write(6,"('Allocated private background t,u,v,z ',
     |  '(lon0:lon1,lat0:lat1,nmonth,nhour)')")

      end subroutine alloc_bgrddata
!-----------------------------------------------------------------------
      end module bgrd_data_module
