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
      use params_module,only:
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
      use input_module,only:  bgrddata_ncfile,
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
        call mkbgrddata(iyear,iday,int(secs),iprint)
      endif
!
      end subroutine get_bgrddata 
! 
!-----------------------------------------------------------------------
      subroutine rd_bgrddata(ncfile,type)
      use input_module,only: mxlen_filename
      implicit none
!
! Args:
      character(len=*),intent(in) :: ncfile,type
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
      write(6,"('Reading GSWM data file ',a)") trim(ncfile)
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
	  ' nmonth=',i4)") nmonth_rd,nmonth
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
	  ' nhour=', i4)") nhour_rd,nhour
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
	  ' nlon=',i4)") nlon_rd,nlon
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
	' nlat=',i4)") nlat_rd,nlat
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
      write(6,"('Completed read from Background data file ',a)") trim(ncfile)
      write(6,"(/,72('-'))")
      end subroutine rd_bgrddata
!-----------------------------------------------------------------------
      end module bgrd_data_module
