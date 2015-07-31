#include "defs.h"
!      
      module amie_module
!      
! Module used to read data from the AMIE outputs (POT,mean energy, 
!   and energy flux).
!
      use params_module,only: nlon,nlat,nlev,nlonp1,nlonp4,nlatp1,
     |  nlevp1,nmlon,nmlonp1,nmlat
      use addfld_module,only: addfld
      implicit none
!
#ifdef SUN
#include "netcdf3.inc"
#else
#include "netcdf.inc"
#endif 
!
! Define parameters for AMIE input data file:
      integer,parameter :: 
     |  mxgdays = 10,     ! maximum number of days of AMIE data
     |  mxtimes = 5881,   ! maximum number of times of AMIE data per day
     |  ithtrns = 30,     ! corresponding to trans lat 40-deg
     |  ithmx = 55,       ! maximum number of latitudes of AMIE data
     |  jmxm = 2*ithmx-1, ! maximum number of global latitudes
     |  lonmx = 36        ! maximum number of longitudes of AMIE data
      integer :: lonp1,latp1
      integer,dimension(mxtimes) :: year,month,day,jday
! Define AMIE output fields
      real :: 
     |  tiepot(nmlonp1,nmlat),tieekv(nmlonp1,nmlat),
     |  tieefx(nmlonp1,nmlat),potg(nlonp1,0:nlatp1),
     |  efxg(nlonp1,0:nlatp1),ekvg(nlonp1,0:nlatp1)
! defined output AMIE fields in TGCM geographic grid
      real,dimension(nlonp4,nlat) :: 
     |  potg_sech, ekvg_sech, efxg_sech
      real,dimension(nmlonp1,-2:nlevp1) :: tiepot_sech
!
! Define fields for AMIE input data file:
! electric potential in Volt
! mean energy in KeV
! energy flux in W/m^2
! amie_cusplat_nh(sh) and amie_cuspmlt_nh(sh) are 
!   AMIE cusp latitude and MLT in NH and SH
! amie_hpi_nh(sh) are AMIE hemi-integrated power
! amie_pcp_nh(sh) are AMIE polar-cap potential drop
! Saved AMIE outputs with suffix _amie
!
      real,allocatable,dimension(:,:,:),save :: ! (lonp1,latp1,ntimes)
     |  amie_pot_nh, amie_pot_sh, amie_ekv_nh, amie_ekv_sh,
     |  amie_efx_nh, amie_efx_sh
      real,allocatable,dimension(:,:),save ::   ! (lonp1,latp1)
     |  pot_nh_amie,pot_sh_amie, ekv_nh_amie,ekv_sh_amie,
     |  efx_nh_amie,efx_sh_amie
      real, allocatable,dimension(:),save :: ! (ntimes)
     |  amie_cusplat_nh, amie_cuspmlt_nh, amie_hpi_nh,
     |  amie_pcp_nh, amie_nh_ut,
     |  amie_cusplat_sh, amie_cuspmlt_sh, amie_hpi_sh,
     |  amie_pcp_sh, amie_sh_ut
      real ::
     |  cusplat_nh_amie, cuspmlt_nh_amie, cusplat_sh_amie,
     |  cuspmlt_sh_amie, hpi_sh_amie, hpi_nh_amie, pcp_sh_amie,
     |  pcp_nh_amie, crad(2), phida(2)
!
      contains 
!-----------------------------------------------------------------------
      subroutine init_amie
!
! Called from tgcm.F 
! (this is not in init.F to avoid circular dependencies)
!
      use input_module,only: amienh,amiesh
!
! Read north hemisphere file:
      if (len_trim(amienh) > 0) then
        write(6,"('Reading AMIENH file ',a)") trim(amienh)
        call rdamie_nh
      endif
!
! Read south hemisphere file:
      if (len_trim(amiesh) > 0) then
        write(6,"('Reading AMIESH file ',a)") trim(amiesh)
        call rdamie_sh
      endif
      end subroutine init_amie
!-----------------------------------------------------------------------
      subroutine rdamie_nh
!
! Read AMIE data for the northern hemisphere from amienh
!
      use input_module,only: amienh
!
! Local:
      character(len=80) :: dskfile
      integer :: istat,ntimes,ncid,ndims,nvars,ngatts,idunlim,ier
      integer :: id_lon,id_lat,id_time,lon,lat,mlt,
     |    idv_year,idv_mon,idv_day,idv_jday,
     |    idv_iter,idv_ut,idv_lon,idv_lat,idv_pot,idv_ekv,
     |    idv_efx,idv_cusplat,idv_cuspmlt,idv_hpi,idv_pcp
!
! Acquire mss file:
      dskfile = ' '
      call getfile(amienh,dskfile)
!
      write(6,"(/,72('-'))")
      write(6,"('RDAMIE_NH: read AMIE data for northern hemisphere:')")
!
! Open netcdf file:
      istat = nf_open(dskfile,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) then
        write(6,"(/,'>>> rdamie_nh: error opening netcdf NH AMIE data ',
     |    'diurnal file ',a)") trim(dskfile)
        call shutdown('rdamie_nh')
      else
        write(6,"('rdamie_nh: opened file ',a)") trim(dskfile)
      endif
!
! Get AMIE grid dimension:
      istat = nf_inq_dimid(ncid,'lon',id_lon)
      istat = nf_inq_dimlen(ncid,id_lon,lonp1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting AMIE longitude dimension')

      istat = nf_inq_dimid(ncid,'lat',id_lat)
      istat = nf_inq_dimlen(ncid,id_lat,latp1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting AMIE latitude dimension')
      write(6,"('lonp1=',i3,' latp1=',i3)") lonp1,latp1
!
! Get time dimension:
      istat = nf_inq_unlimdim(ncid,id_time)
      istat = nf_inq_dimlen(ncid,id_time,ntimes)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting time dimension')
!     write(6,"('rdamie_nh: ntimes=',i3)") ntimes
! 
! Search for requested AMIE output fields
      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
!
! Get 1-D AMIE fields (ntimes)
      istat = nf_inq_varid(ncid,'year',idv_year)
      istat = nf_get_var_int(ncid,idv_year,year)
      istat = nf_inq_varid(ncid,'month',idv_mon)
      istat = nf_get_var_int(ncid,idv_mon,month)
      istat = nf_inq_varid(ncid,'day',idv_day)
      istat = nf_get_var_int(ncid,idv_day,day)
      istat = nf_inq_varid(ncid,'jday',idv_jday)
      istat = nf_get_var_int(ncid,idv_jday,jday)
!
! Allocate 1-d fields:
      if (.not. allocated(amie_nh_ut))
     |  allocate(amie_nh_ut(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_nh_ut: ntimes=',i3)")ntimes
      if (.not. allocated(amie_cusplat_nh))
     |  allocate(amie_cusplat_nh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_cusplat_nh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_cuspmlt_nh))
     |  allocate(amie_cuspmlt_nh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_cuspmlt_nh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_hpi_nh))
     |  allocate(amie_hpi_nh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_hpi_nh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_pcp_nh))
     |  allocate(amie_pcp_nh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_pcp_nh: ntimes=',i3)")ntimes
!
! Get ut
      istat = nf_inq_varid(ncid,'ut',idv_ut)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE UT id')
      istat = nf_get_var_double(ncid,idv_ut,amie_nh_ut)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable ut')
!
! Get HPI
      istat = nf_inq_varid(ncid,'hpi',idv_hpi)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE hpi id')
      istat = nf_get_var_double(ncid,idv_hpi,amie_hpi_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable hpi')
!
! Get PCP
      istat = nf_inq_varid(ncid,'pcp',idv_pcp)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE pcp id')
      istat = nf_get_var_double(ncid,idv_pcp,amie_pcp_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable pcp')
!
! Get cusplat
      istat = nf_inq_varid(ncid,'cusplat',idv_cusplat)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE cusplat id')
      istat = nf_get_var_double(ncid,idv_cusplat,amie_cusplat_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable cusplat')
!
! Get cuspmlt
      istat = nf_inq_varid(ncid,'cuspmlt',idv_cuspmlt)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE cusplat id')
      istat = nf_get_var_double(ncid,idv_cuspmlt,amie_cuspmlt_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable cuspmlt')
!
! Allocate 2-d fields:
      if (.not. allocated(pot_nh_amie))
     |   allocate(pot_nh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' pot_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
      if (.not. allocated(ekv_nh_amie))
     |   allocate(ekv_nh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' ekv_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
      if (.not. allocated(efx_nh_amie))
     |   allocate(efx_nh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' efx_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
!
! Allocate 3-d fields:
      if (.not. allocated(amie_pot_nh)) 
     |   allocate(amie_pot_nh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_pot_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
      if (.not. allocated(amie_ekv_nh)) 
     |   allocate(amie_ekv_nh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_ekv_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
      if (.not. allocated(amie_efx_nh)) 
     |   allocate(amie_efx_nh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_nh: error allocating',
     |  ' amie_efx_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
!
! Get 3-D AMIE fields (lon,lat,ntimes)
!
! AMIE electric potential
      istat = nf_inq_varid(ncid,'pot',idv_pot)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE electric potential id')
      istat = nf_get_var_double(ncid,idv_pot,amie_pot_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable pot')
!
! AMIE mean energy
      istat = nf_inq_varid(ncid,'ekv',idv_ekv)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE mean energy id')
      istat = nf_get_var_double(ncid,idv_ekv,amie_ekv_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable ekv')
!
! AMIE energy flux
      istat = nf_inq_varid(ncid,'efx',idv_efx)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE energy flux id')
      istat = nf_get_var_double(ncid,idv_efx,amie_efx_nh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_nh: Error getting NH AMIE variable efx')
!
! Close the file:
      istat = nf_close(ncid)
      write(6,"('Completed read from NH AMIE data file ',a)")
     |   trim(dskfile)
      write(6,"(72('-'),/)")
      end subroutine rdamie_nh
!-----------------------------------------------------------------------
      subroutine rdamie_sh
!
! Read AMIE data for the southern hemisphere from amiesh
!
      use input_module,only: amiesh
!
! Local:
      character(len=80) :: dskfile
      integer :: istat,ntimes,ncid,ndims,nvars,ngatts,idunlim,ier
      integer :: id_lon,id_lat,id_time,lon,lat,mlt,
     |    idv_year,idv_mon,idv_day,idv_jday,
     |    idv_iter,idv_ut,idv_lon,idv_lat,idv_pot,idv_ekv,
     |    idv_efx,idv_cusplat,idv_cuspmlt,idv_hpi,idv_pcp
!
! Acquire mss file:
      dskfile = ' '
      call getfile(amiesh,dskfile)
!
      write(6,"(/,72('-'))")
      write(6,"('RDAMIE_SH: read AMIE data for southern hemisphere:')")
!
! Open netcdf file:
      istat = nf_open(dskfile,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) then
        write(6,"(/,'>>> rdamie_nh: error opening netcdf SH AMIE data ',
     |    'diurnal file ',a)") trim(dskfile)
        call shutdown('rdamie_sh')
      endif
!
! Get AMIE grid dimension:
      istat = nf_inq_dimid(ncid,'lon',id_lon)
      istat = nf_inq_dimlen(ncid,id_lon,lonp1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting AMIE longitude dimension')
      istat = nf_inq_dimid(ncid,'lat',id_lat)
      istat = nf_inq_dimlen(ncid,id_lat,latp1)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting AMIE latitude dimension')
!
! Get time dimension:
      istat = nf_inq_unlimdim(ncid,id_time)
      istat = nf_inq_dimlen(ncid,id_time,ntimes)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting time dimension')
! 
! Search for requested AMIE output fields
      istat = nf_inq(ncid,ndims,nvars,ngatts,idunlim)
!
! Get 1-D AMIE fields (ntimes)
      istat = nf_inq_varid(ncid,'year',idv_year)
      istat = nf_get_var_int(ncid,idv_year,year)
      istat = nf_inq_varid(ncid,'month',idv_mon)
      istat = nf_get_var_int(ncid,idv_mon,month)
      istat = nf_inq_varid(ncid,'day',idv_day)
      istat = nf_get_var_int(ncid,idv_day,day)
      istat = nf_inq_varid(ncid,'jday',idv_jday)
      istat = nf_get_var_int(ncid,idv_jday,jday)
!
! Allocate 1-d fields:
      if (.not. allocated(amie_sh_ut))
     |   allocate(amie_sh_ut(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_sh_ut: ntimes=',i3)")ntimes
      if (.not. allocated(amie_cusplat_sh))
     |   allocate(amie_cusplat_sh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_cusplat_sh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_cuspmlt_sh))
     |   allocate(amie_cuspmlt_sh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_cuspmlt_sh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_hpi_sh))
     |   allocate(amie_hpi_sh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_hpi_sh: ntimes=',i3)")ntimes
      if (.not. allocated(amie_pcp_sh))
     |   allocate(amie_pcp_sh(ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_pcp_sh: ntimes=',i3)")ntimes
!
! Get ut
      istat = nf_inq_varid(ncid,'ut',idv_ut)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE UT id')
      istat = nf_get_var_double(ncid,idv_ut,amie_sh_ut)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable ut')
!
! Get HPI
      istat = nf_inq_varid(ncid,'hpi',idv_hpi)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE hpi id')
      istat = nf_get_var_double(ncid,idv_hpi,amie_hpi_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable hpi')
!
! Get PCP
      istat = nf_inq_varid(ncid,'pcp',idv_pcp)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE pcp id')
      istat = nf_get_var_double(ncid,idv_pcp,amie_pcp_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable hpi')
!
! Get cusplat
      istat = nf_inq_varid(ncid,'cusplat',idv_cusplat)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE cusplat id')
      istat = nf_get_var_double(ncid,idv_cusplat,amie_cusplat_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable cusplat')
!
! Get cuspmlt
      istat = nf_inq_varid(ncid,'cuspmlt',idv_cuspmlt)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE cusplat id')
      istat = nf_get_var_double(ncid,idv_cuspmlt,amie_cuspmlt_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable cuspmlt')
!
! Allocate 2-d fields:
      if (.not. allocated(pot_sh_amie))
     |   allocate(pot_sh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' pot_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
      if (.not. allocated(ekv_sh_amie))
     |   allocate(ekv_sh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' ekv_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
      if (.not. allocated(efx_sh_amie))
     |   allocate(efx_sh_amie(lonp1,latp1),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' efx_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1

!
! Allocate 3-d fields:
      if (.not. allocated(amie_pot_sh))
     |   allocate(amie_pot_sh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_pot_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
      if (.not. allocated(amie_ekv_sh))
     |   allocate(amie_ekv_sh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_ekv_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
      if (.not. allocated(amie_efx_sh))
     |   allocate(amie_efx_sh(lonp1,latp1,ntimes),stat=ier)
      if (ier /= 0) write(6,"('>>> rdamie_sh: error allocating',
     |  ' amie_efx_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)")
     |  lonp1,latp1,ntimes
!
! Get 3-D AMIE fields (lon,lat,ntimes)
!
! AMIE electric potential
      istat = nf_inq_varid(ncid,'pot',idv_pot)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE electric potential id')
      istat = nf_get_var_double(ncid,idv_pot,amie_pot_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable pot')
!
! AMIE mean energy
      istat = nf_inq_varid(ncid,'ekv',idv_ekv)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE mean energy id')
      istat = nf_get_var_double(ncid,idv_ekv,amie_ekv_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable ekv')
!
! AMIE energy flux
      istat = nf_inq_varid(ncid,'efx',idv_efx)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE energy flux id')
      istat = nf_get_var_double(ncid,idv_efx,amie_efx_sh)
      if (istat /= NF_NOERR) call rpt_ncerr(istat,
     |  'rdamie_sh: Error getting SH AMIE variable efx')
!
! Close the file:
      istat = nf_close(ncid)
      write(6,"('Completed read from SH AMIE data file ',a)")
     |   trim(dskfile)
      write(6,"(72('-'),/)")
      end subroutine rdamie_sh
!-----------------------------------------------------------------------
      subroutine getamie(iyear,iday,iutsec,amie_ibkg,iprint)
!
! Read AMIE outputs from amie_ncfile file, returning electric potential,
! auroral mean energy and energy flux at current date and time,
! and the data is linearly interpolated to the model time
! gl - 12/07/2002
!
      use hist_module,only: modeltime
      use input_module,only: iamie,amiesh,amienh
      use cons_module,only: pi,ylonm,ylatm
      use magfield_module,only: im,jm,dim,djm
!
! Args:
      integer,intent(in) :: iyear,iday,iutsec,amie_ibkg,iprint
!
! Local:
      real :: 
     |   potm(lonp1,jmxm),efxm(lonp1,jmxm),ekvm(lonp1,jmxm),
     |   alat(jmxm),alon(lonp1),alatm(jmxm),alonm(lonp1)
      integer :: ier,lw,liw,isign,intpol(2)
      integer,allocatable :: iw(:)
      real,allocatable :: w(:)
      integer :: i,j,ndmon_nl(13),ndmon_lp(13),nonleap
      integer :: nn, iset, iset1, m, mp1, n, amie_debug,amieday
      integer :: iboxcar,inpt
      real    :: model_ut, denoma, f1, f2, fmin, fmax, maxefx
      real    :: del,xmlt,dmlat,dlatm,dlonm,dmltm,rot,dtr
!   
! Data: GSWM data given at the middle of each month 
!       -> day of year for this data assuming non-leap year
!         J  F  M  A   M   J   J   A   S   O   N   D   J
      data ndmon_nl                                ! non_leap
     |  /31,59,90,120,151,181,212,243,273,304,334,365,396/ 
      data ndmon_lp                                ! leap year
     |  /31,60,91,121,152,182,213,244,274,305,335,366,397/ 
!
! Reading Charley Barth's auroral modification factor
      integer, parameter:: nxtimes=47
      integer :: nxyear,nxday(nxtimes)
      real :: fxfact_nh, fxfact_sh,fxtime,fxtime1,fxtime2
      real, dimension(nxtimes):: ut_nh,ut_sh,fx_nh,fx_sh
! 
! External:
      real,external :: finterp
!
      nxyear = 1998
      nxday(1:nxtimes)=(/267,
     |  267, 267, 267, 267, 267, 267, 267,
     |  267, 267, 267, 267, 267, 267, 267,
     |  267, 268, 268, 268, 268, 268, 268, 
     |  268, 268, 268, 268, 268, 268, 268, 
     |  268, 268, 269, 269, 269, 269, 269, 
     |  269, 269, 269, 269, 269, 269, 269, 
     |  269, 269, 269, 270/)
      ut_nh(1:nxtimes)=(/0.0,
     |   0.68872, 2.28569, 3.88347, 5.47933, 7.07654, 8.67133, 
     |   10.2672, 11.8640, 13.4619, 15.0575, 16.6542, 18.2520, 
     |   19.8473, 21.4437, 23.0409, 0.63913, 2.23508, 3.83221, 
     |   5.42711, 7.02325, 8.62048, 10.2171, 11.8147, 13.4118, 
     |   15.0064, 16.6022, 18.1989, 19.7967, 21.3921, 22.9884, 
     |   0.58565, 2.18393, 3.78000, 5.37563, 6.97244, 8.56865, 
     |   10.1659, 11.7654, 13.3620, 14.9582, 16.5557, 18.1509,  
     |   19.7470, 21.3442, 22.9422,  0./)
      fx_nh(1:nxtimes)=(/1.0,
     |  0.954679, 1.01677, 0.805276, 0.662938, 0.680216, 1.09095, 
     |  0.883288,0.972496, 0.720709, 0.579112, 0.818075, 1.28912, 
     |   1.35952, 1.36740,  1.03170,  1.34414,  1.53157, 0.90830, 
     |  0.824193, 1.10184,  1.35307, 0.999406, 0.788716, 0.76273, 
     |  0.647832, 0.57241,  0.52249,  1.48004,  2.36089, 1.86270,  
     |   1.87322, 1.80098,  1.57992,  1.33901,  1.21373, 1.42878, 
     |   1.28782,  1.2731,  1.31435,  1.13873, 0.984014, 1.04712,
     |   1.48459, 1.75276,  1.94827,  1./)

      ut_sh(1:nxtimes)=(/0.0,
     |  0.104409, 1.70090,  3.29847,  4.89378,  6.49047, 8.08820,  
     |   9.68356, 11.2801,  12.8777,  14.4745,  16.0726, 17.6666, 
     |   19.2649, 20.8609,  22.4578,  0.05568,  1.65123, 3.24788,  
     |   4.84236, 6.44144,  8.03825,  9.63280,  11.2301, 12.8267,  
     |   14.4227, 16.0197,  17.6161,  19.2119,  20.8103, 22.4062,  
     |   0.02522, 1.59759,  3.19667,  4.79183,  6.38823, 7.98401, 
     |   9.58092, 11.1789,  12.7793,  14.3752,  15.9723, 17.5671,
     |   19.1629, 20.7597,  22.3572,  0./)

      fx_sh(1:nxtimes)=(/1.0,
     |  0.544252, 0.487370, 0.572375, 0.605849, 0.605344, 0.776822, 
     |  0.895618, 0.641604, 0.590224, 0.638560, 0.739044, 0.863948,  
     |   1.22977,  1.16256, 0.881014, 0.669191, 0.617876, 0.498846, 
     |  0.998803, 0.826644, 0.899562,  1.41258,  1.05342, 0.972168,  
     |   1.00838,  1.50718,  2.00560,  2.01770,  1.85613,  1.31682, 
     |  0.920968,  1.64630,  2.33071,  1.96759,  1.83485,  1.79565, 
     |   1.74479,  1.52049,  1.76034,  1.85011,  1.32154,  1.19636,
     |   1.28152, 1.42212,  1.49638,  1./)
!
      if (iprint > 0) then
        write(6,"(/,72('-'))")
        write(6,"('GETAMIE:')")
        write(6,"('Modeltime(1:4)= ',4i5)") modeltime(1:4)
        write(6,"('Initial requested iyear=',i4,' iday=',i3,' iutsec=',
     |    i10)") iyear,iday,iutsec
      endif
! 
! Check for leap year
      nonleap = int(mod(real(iyear),4.))  ! nonleap year /= 0
!     if (iyear.eq.2000) nonleap = 1
! 
! Check times:
!
      if (year(1) /= iyear) then
        write(6,"('>>> getamie: wrong year,iyear=',2i4)")
     |     year(1),iyear
        call shutdown('getamie')
      endif
      amieday = jday(1)
      if (amie_ibkg > 0) amieday = jday(1) - 1
!     if (amieday /= iday .and. amie_ibkg == 0) then
!       write(6,"('>>> getamie: wrong amieday,modelday=',2i4)")
!    |     amieday,iday
!       call shutdown('getamie')
!     endif
      if (modeltime(3).lt.0.or.modeltime(3).gt.59) then
        write(6,"('>>> getamie: bad imin=',3i4)") modeltime(3)
        call shutdown('getamie')
      endif
      if (modeltime(4).lt.0.or.modeltime(4).gt.59) then
        write(6,"('>>> getamie: bad isec=',3i4)") modeltime(4)
        call shutdown('getamie')
      endif
!     model_ut = (modeltime(1)-amieday)*24. + float(iutsec)/3600.
      model_ut = float(iutsec)/3600.

      maxefx=50.
      if (iyear==2004 .and. amieday==314 .and. abs(model_ut-7.333)<0.1) 
     |      maxefx=21.
      if (iyear==2004 .and. amieday==314 .and. abs(model_ut-8.833)<0.1) 
     |      maxefx=14.
      if (iyear==2004 .and. amieday==318) maxefx=10.
      if (iyear==2004 .and. amieday==318 .and. model_ut < 5.75) 
     |      maxefx=3.
      if (iyear==2004 .and. amieday==318 .and. abs(model_ut-4.75)<0.1) 
     |      maxefx=1.
      if (iyear==2004 .and. amieday==318 .and. abs(model_ut-21.0)<0.1) 
     |      maxefx=2.

! Find out the multiplication factor from Charley Barth
      fxfact_nh = 1.
      fxfact_sh = 1.
      if (nxyear == iyear) then
! For the southern hemisphere
	iset = 0
	iset1 = nxtimes
        do i=1,nxtimes
          fxtime=(nxday(i)-iday)*24.+ut_nh(i)
	  if (fxtime < model_ut) iset=i
	enddo
	iset1 = iset + 1
	if (iset==0) iset=1
	if (iset==nxtimes) iset1=iset
        fxtime1=(iday-nxday(iset))*24.+ut_nh(iset)
        fxtime2=(iday-nxday(iset1))*24.+ut_nh(iset1)
	denoma = fxtime2 - fxtime1
	if (denoma == 0.) then
	  f1 = 1.
	  f2 = 0.
	else
!	  f1 = (fxtime2-model_ut)/denoma
!	  f2 = (model_ut-fxtime1)/denoma
	  f1 = (fxtime2-(model_ut+(modeltime(1)-jday(i))*24.))/denoma
	  f2 = (model_ut+(modeltime(1)-jday(i))*24.-fxtime1)/denoma
	endif
	fxfact_nh = f1*fx_nh(iset)+f2*fx_nh(iset1)
! For the southern hemisphere
	iset = nxtimes
	iset1 = nxtimes
        do i=1,nxtimes
          fxtime=(nxday(i)-iday)*24.+ut_sh(i)
!	  if (fxtime < model_ut) iset=i
	  if (fxtime < model_ut+(modeltime(1)-jday(i))*24.) iset=i
	enddo
	iset1 = iset + 1
	if (iset==nxtimes) iset1=iset
        fxtime1=(iday-nxday(iset))*24.+ut_sh(iset)
        fxtime2=(iday-nxday(iset1))*24.+ut_sh(iset1)
	denoma = fxtime2 - fxtime1
	if (denoma == 0.) then
	  f1 = 1.
	  f2 = 0.
	else
!	  f1 = (fxtime2-model_ut)/denoma
!	  f2 = (model_ut-fxtime1)/denoma
	  f1 = (fxtime2-(model_ut+(modeltime(1)-jday(i))*24.))/denoma
	  f2 = (model_ut+(modeltime(1)-jday(i))*24.-fxtime1)/denoma
	endif
	fxfact_sh = f1*fx_sh(iset)+f2*fx_sh(iset1)
      endif
! 
! Interpolate AMIE data to modeltime iutsec
! amie_ibkg = 0  use real UT AMIE data
!      = 1  use the first AMIE volumne as the background
!      = 2  use the 24-hr average AMIE volumne as the background
      pot_sh_amie(:,:) = 0.
      ekv_sh_amie(:,:) = 0.
      efx_sh_amie(:,:) = 0.
      cusplat_sh_amie = 0.
      cuspmlt_sh_amie = 0.
      hpi_sh_amie = 0.
      pcp_sh_amie = 0.
!
      n = len_trim(amiesh)

      iboxcar = 0
      if (n > 0) then
        nn = size(amie_sh_ut)
        if (jday(nn) < modeltime(1)) then
          write(6,"('getamie: Model date beyond the AMIE last Data: ',
     |   'modelday, amieday =',2i5)") modeltime(1),jday(nn)
          iamie = 0
          return
        endif
        if (amie_ibkg == 0) then
          iset = nn
          iset1 = nn
          do i=1,nn
!           if (amie_sh_ut(i) < model_ut) iset = i
            if (amie_sh_ut(i) < model_ut+(modeltime(1)-jday(i))*24.) 
     |          iset = i
          enddo
!       write(6,"('getamie: AMIE SH Data nn,iset,day1,day2=',4i5)")
!    |    nn,iset,jday(1),jday(nn)
          iset1 = iset + 1
          if (iset == nn) iset1 = iset
! get rid of outragous flux values 
      inpt = 0
      do j=1,latp1
        do i=1,lonp1 
          if (amie_efx_sh(i,j,iset) > maxefx) then
             if (i > 1) then
               amie_efx_sh(i,j,iset)=amie_efx_sh(i-1,j,iset)
             else 
               amie_efx_sh(i,j,iset)=amie_efx_sh(i+1,j,iset)
             endif
             if (j < 12) amie_efx_sh(i,j,iset) = 1.
             inpt = inpt + 1
          endif
          if (iset1 /= iset .and. amie_efx_sh(i,j,iset1) > maxefx) then
             if (i > 1) then
               amie_efx_sh(i,j,iset1)=amie_efx_sh(i-1,j,iset1)
             else 
               amie_efx_sh(i,j,iset1)=amie_efx_sh(i+1,j,iset1)
             endif
             if (j < 12) amie_efx_sh(i,j,iset1) = 1.
             inpt = inpt + 1
          endif
        enddo
      enddo
      if (inpt > 0) write(6,"('maxefx = ',f4.1,
     |	 '  Bad efx points SH= ',i3)")maxefx,inpt

          denoma = amie_sh_ut(iset1) - amie_sh_ut(iset)
          if (denoma == 0.) then
            f1 = 1.
            f2 = 0.
          else
!           f1 = (amie_sh_ut(iset1) - model_ut)/denoma
!           f2 = (model_ut - amie_sh_ut(iset))/denoma
            f1 = (amie_sh_ut(iset1) - (model_ut+(modeltime(1)-
     |             jday(iset1))*24.))/denoma
            f2 = (model_ut+(modeltime(1)-jday(iset1))*24. - 
     |             amie_sh_ut(iset))/denoma
          endif
!       write(6,"('getamie: AMIE SH Data n,iset,modeltime,f1,f2 =',
!    |    4i5,2f5.2)")n,iset,modeltime(1),jday(iset1),f1,f2
!       write(6,"('getamie: AMIE SH Data model_day,model_ut,amie_day,',
!    |    'amie_ut,f1,f2,iset,iset1 =',i4,f7.1,i4,f7.1,2f5.2,2i3)")
!    |    modeltime(1),model_ut,jday(iset),amie_sh_ut(iset),f1,f2,
!    |    iset,iset1
          cusplat_sh_amie = (f1*amie_cusplat_sh(iset1) + 
     |                        f2*amie_cusplat_sh(iset))
          cuspmlt_sh_amie = (f1*amie_cuspmlt_sh(iset1) + 
     |                        f2*amie_cuspmlt_sh(iset))
          hpi_sh_amie = (f1*amie_hpi_sh(iset1) + f2*amie_hpi_sh(iset))
          pcp_sh_amie = (f1*amie_pcp_sh(iset1) + f2*amie_pcp_sh(iset))
          if (iboxcar == 0) then
            pot_sh_amie(:,:) = (f1*amie_pot_sh(:,:,iset1) + 
     |                        f2*amie_pot_sh(:,:,iset))
            ekv_sh_amie(:,:) = (f1*amie_ekv_sh(:,:,iset1) + 
     |                        f2*amie_ekv_sh(:,:,iset))
            efx_sh_amie(:,:) = (f1*amie_efx_sh(:,:,iset1) + 
     |                        f2*amie_efx_sh(:,:,iset))
          else
            call boxcar_ave(amie_pot_sh,pot_sh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
            call boxcar_ave(amie_efx_sh,efx_sh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
            call boxcar_ave(amie_ekv_sh,ekv_sh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
          endif
        else 
          if (amie_ibkg == 1) then
             pot_sh_amie(:,:) = amie_pot_sh(:,:,1)
             ekv_sh_amie(:,:) = amie_ekv_sh(:,:,1)
             efx_sh_amie(:,:) = amie_efx_sh(:,:,1)
             cusplat_sh_amie = amie_cusplat_sh(1)
             cuspmlt_sh_amie = amie_cuspmlt_sh(1)
             hpi_sh_amie = amie_hpi_sh(1)
             pcp_sh_amie = amie_pcp_sh(1)
          else if (amie_ibkg == 3) then
             pot_sh_amie(:,:) = amie_pot_sh(:,:,241)
             ekv_sh_amie(:,:) = amie_ekv_sh(:,:,241)
             efx_sh_amie(:,:) = amie_efx_sh(:,:,241)
             cusplat_sh_amie = amie_cusplat_sh(241)
             cuspmlt_sh_amie = amie_cuspmlt_sh(241)
             hpi_sh_amie = amie_hpi_sh(241)
             pcp_sh_amie = amie_pcp_sh(241)
           else
             do i=1,nn
               pot_sh_amie(:,:) = pot_sh_amie(:,:) + amie_pot_sh(:,:,1)
               ekv_sh_amie(:,:) = ekv_sh_amie(:,:) + amie_ekv_sh(:,:,1)
               efx_sh_amie(:,:) = efx_sh_amie(:,:) + amie_efx_sh(:,:,1)
               cusplat_sh_amie = cusplat_sh_amie + amie_cusplat_sh(1)
               cuspmlt_sh_amie = cuspmlt_sh_amie + amie_cuspmlt_sh(1)
               hpi_sh_amie = hpi_sh_amie + amie_hpi_sh(1)
               pcp_sh_amie = pcp_sh_amie + amie_pcp_sh(1)
             enddo
             pot_sh_amie(:,:) = pot_sh_amie(:,:)/nn
             ekv_sh_amie(:,:) = ekv_sh_amie(:,:)/nn
             efx_sh_amie(:,:) = efx_sh_amie(:,:)/nn
             cusplat_sh_amie = cusplat_sh_amie/nn
             cuspmlt_sh_amie = cuspmlt_sh_amie/nn
             hpi_sh_amie = hpi_sh_amie/nn
             pcp_sh_amie = pcp_sh_amie/nn
           endif
         endif
      endif
! get rid of outragous flux values 
      inpt = 0
      do j=1,latp1
        do i=1,lonp1 
          if (efx_sh_amie(i,j) > maxefx) then
!            efx_sh_amie(i,j) = 1.
             if (i > 1) then
               efx_sh_amie(i,j)=efx_sh_amie(i-1,j)
             else 
               efx_sh_amie(i,j)=efx_sh_amie(i+1,j)
             endif
             inpt = inpt + 1
          endif
        enddo
      enddo
      if (inpt > 0) write(6,"('Total bad efx points SH= ',i3)")inpt
! 
! get NH AMIE data
      pot_nh_amie(:,:) = 0.
      ekv_nh_amie(:,:) = 0.
      efx_nh_amie(:,:) = 0.
      cusplat_nh_amie = 0.
      cuspmlt_nh_amie = 0.
      hpi_nh_amie = 0.
      pcp_nh_amie = 0.

      n = len_trim(amienh)
      iboxcar = 0
      if (n > 0) then
        nn = size(amie_nh_ut)
!       if (amie_nh_ut(nn) < model_ut+(iday-jday(n))*24.) then
        if (jday(nn) < modeltime(1)) then
          write(6,"('getamie: Model date beyond the AMIE last Data: ',
     |   'modelday, amieday =',2i5)") modeltime(1),jday(nn)
          iamie = 0
          return
        endif
!       write(6,"('getamie: Interpolate AMIE NH Data nn=',i3)")nn
        if (amie_ibkg == 0) then
          iset = 0
          iset1 = nn
          do i=1,nn
            if (amie_nh_ut(i) < model_ut+(modeltime(1)-jday(i))*24.) 
     |          iset = i
          enddo
          iset1 = iset + 1
          if (iset == 0) iset = 1
          if (iset == nn) iset1 = iset
! get rid of outragous flux values 
      inpt = 0
      do j=1,latp1
        do i=1,lonp1 
          if (amie_efx_nh(i,j,iset) > maxefx) then
             if (i > 1) then
               amie_efx_nh(i,j,iset)=amie_efx_nh(i-1,j,iset)
             else 
               amie_efx_nh(i,j,iset)=amie_efx_nh(i+1,j,iset)
             endif
             if (j < 12) amie_efx_nh(i,j,iset) = 1.
             inpt = inpt + 1
          endif
          if (iset1 /= iset .and. amie_efx_nh(i,j,iset1) > maxefx) then
             if (i > 1) then
               amie_efx_nh(i,j,iset1)=amie_efx_nh(i-1,j,iset1)
             else 
               amie_efx_nh(i,j,iset1)=amie_efx_nh(i+1,j,iset1)
             endif
             if (j < 12) amie_efx_nh(i,j,iset1) = 1.
             inpt = inpt + 1
          endif
        enddo
      enddo
      if (inpt > 0) write(6,"('maxefx = ',f4.1,
     |	 '  Bad efx points NH= ',i3)")maxefx,inpt

          denoma = amie_nh_ut(iset1) - amie_nh_ut(iset)
          if (denoma == 0.) then
            f1 = 1.
            f2 = 0.
          else
!           f1 = (amie_nh_ut(iset1) - model_ut)/denoma
!           f2 = (model_ut - amie_nh_ut(iset))/denoma
            f1 = (amie_nh_ut(iset1) - (model_ut+(modeltime(1)-
     |             jday(iset1))*24.))/denoma
            f2 = (model_ut+(modeltime(1)-jday(iset1))*24. - 
     |             amie_nh_ut(iset))/denoma
          endif
!       write(6,"('getamie: AMIE NH Data model_day,model_ut,amie_day,',
!    |    'amie_ut,f1,f2,iset,iset1 =',i4,f7.1,i4,f7.1,2f5.2,2i3)")
!    |    modeltime(1),model_ut,jday(iset),amie_nh_ut(iset),f1,f2,
!    |    iset,iset1
!
          cusplat_nh_amie = (f1*amie_cusplat_nh(iset1) + 
     |                      f2*amie_cusplat_nh(iset))
          cuspmlt_nh_amie = (f1*amie_cuspmlt_nh(iset1) + 
     |                      f2*amie_cuspmlt_nh(iset))
          hpi_nh_amie = (f1*amie_hpi_nh(iset1) + f2*amie_hpi_nh(iset))
          pcp_nh_amie = (f1*amie_pcp_nh(iset1) + f2*amie_pcp_nh(iset))
          if (iboxcar == 0) then
            pot_nh_amie(:,:) = (f1*amie_pot_nh(:,:,iset1) + 
     |                        f2*amie_pot_nh(:,:,iset))
            ekv_nh_amie(:,:) = (f1*amie_ekv_nh(:,:,iset1) + 
     |                        f2*amie_ekv_nh(:,:,iset))
            efx_nh_amie(:,:) = (f1*amie_efx_nh(:,:,iset1) + 
     |                        f2*amie_efx_nh(:,:,iset))
!          write(6,"('ekv_nh_amie min, max = ',2e12.4)")
!    |       minval(ekv_nh_amie),maxval(ekv_nh_amie)
          else
            call boxcar_ave(amie_pot_nh,pot_nh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
!           call fminmax(amie_pot_nh(:,:,iset),lonp1*latp1,fmin,fmax)
!           write(6,"('AMIE pot max,min = ',2f8.0)")fmax,fmin
!           call fminmax(pot_nh_amie(:,:),lonp1*latp1,fmin,fmax)
!           write(6,"('boxcar_ave AMIE pot max,min= ',2f8.0)")fmax,fmin
            call boxcar_ave(amie_efx_nh,efx_nh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
!           call fminmax(amie_efx_nh(:,:,iset),lonp1*latp1,fmin,fmax)
!           write(6,"('AMIE efx max,min = ',2f8.0)")fmax,fmin
!           call fminmax(efx_nh_amie(:,:),lonp1*latp1,fmin,fmax)
!           write(6,"('boxcar_ave AMIE efx max,min= ',2f8.0)")fmax,fmin
            call boxcar_ave(amie_ekv_nh,ekv_nh_amie,lonp1,latp1,
     |             nn,iset,iboxcar)
!           call fminmax(amie_ekv_nh(:,:,iset),lonp1*latp1,fmin,fmax)
!           write(6,"('AMIE ekv max,min = ',2f8.0)")fmax,fmin
!           call fminmax(ekv_nh_amie(:,:),lonp1*latp1,fmin,fmax)
!           write(6,"('boxcar_ave AMIE ekv max,min= ',2f8.0)")fmax,fmin
          endif
        else 
          if (amie_ibkg == 1) then
             pot_nh_amie(:,:) = amie_pot_nh(:,:,1)
             ekv_nh_amie(:,:) = amie_ekv_nh(:,:,1)
             efx_nh_amie(:,:) = amie_efx_nh(:,:,1)
             cusplat_nh_amie = amie_cusplat_nh(1)
             cuspmlt_nh_amie = amie_cuspmlt_nh(1)
             hpi_nh_amie = amie_hpi_nh(1)
             pcp_nh_amie = amie_pcp_nh(1)
          else if (amie_ibkg == 3) then
             pot_nh_amie(:,:) = amie_pot_nh(:,:,241)
             ekv_nh_amie(:,:) = amie_ekv_nh(:,:,241)
             efx_nh_amie(:,:) = amie_efx_nh(:,:,241)
             cusplat_nh_amie = amie_cusplat_nh(241)
             cuspmlt_nh_amie = amie_cuspmlt_nh(241)
             hpi_nh_amie = amie_hpi_nh(241)
             pcp_nh_amie = amie_pcp_nh(241)
           else
             do i=1,nn
               pot_nh_amie(:,:) = pot_nh_amie(:,:) + amie_pot_nh(:,:,1)
               ekv_nh_amie(:,:) = ekv_nh_amie(:,:) + amie_ekv_nh(:,:,1)
               efx_nh_amie(:,:) = efx_nh_amie(:,:) + amie_efx_nh(:,:,1)
               cusplat_nh_amie = cusplat_nh_amie + amie_cusplat_nh(1)
               cuspmlt_nh_amie = cuspmlt_nh_amie + amie_cuspmlt_nh(1)
               hpi_nh_amie = hpi_nh_amie + amie_hpi_nh(1)
               pcp_nh_amie = pcp_nh_amie + amie_pcp_nh(1)
             enddo
             pot_nh_amie(:,:) = pot_nh_amie(:,:)/nn
             ekv_nh_amie(:,:) = ekv_nh_amie(:,:)/nn
             efx_nh_amie(:,:) = efx_nh_amie(:,:)/nn
             cusplat_nh_amie = cusplat_nh_amie/nn
             cuspmlt_nh_amie = cuspmlt_nh_amie/nn
             hpi_nh_amie = hpi_nh_amie/nn
             pcp_nh_amie = pcp_nh_amie/nn
           endif
         endif
      endif
! get rid of outragous flux values 
      inpt = 0
      do j=1,latp1
        do i=1,lonp1 
          if (efx_nh_amie(i,j) > maxefx) then
!            efx_nh_amie(i,j) = 1.
             if (i > 1) then
               efx_nh_amie(i,j)=efx_nh_amie(i-1,j)
             else 
               efx_nh_amie(i,j)=efx_nh_amie(i+1,j)
             endif
             inpt = inpt + 1
          endif
! Reduce auroral flux intensity for HPI < 80 GW for Sep. 2005 event
!         if (hpi_nh_amie < 70.) efx_nh_amie(i,j) = efx_nh_amie(i,j)/2.
        enddo
      enddo
      if (inpt > 0) write(6,"('Total bad efx points NH = ',i3)")inpt
! Multiply Barth's factor for auroral fluxes
      if (nxyear==iyear) then
        do j=1,19
          do i=1,lonp1
            if (ekv_nh_amie(i,j) > 2.5)
     |               ekv_nh_amie(i,j) = 1.5*ekv_nh_amie(i,j)
            if (ekv_sh_amie(i,j) > 2.5)
     |               ekv_sh_amie(i,j) = 1.5*ekv_sh_amie(i,j)
          enddo
        enddo
      endif
      fxfact_nh = 1.
      fxfact_sh = 1.
      if (fxfact_nh > 1. .or. fxfact_sh > 1.) then
        write(6,"('Barth Multiplication factors = ',2i5,2f7.3)")
     |    iyear, nxyear, fxfact_nh, fxfact_sh
        efx_nh_amie(:,:) = efx_nh_amie(:,:)*fxfact_nh
        efx_sh_amie(:,:) = efx_sh_amie(:,:)*fxfact_sh
      endif
!
!  The OLTMAX latitude also defines the co-latitude theta0, which in
!   turn determines crit1(+2.5deg) and crit2(-12.5deg) which are used
!   in TIE-GCM as the boundaries of the polar cap and the region of
!   influence of the high-lat potential versus the low-lat dynamo potential
!  Define this latitude to be between 70 and 77.5 degrees
!
!     if (cusplat_sh_amie > 65.0) then
!       cusplat_sh_amie = 65.0
!       cuspmlt_sh_amie = 11.
!     endif
      if (cusplat_sh_amie > 75.0) then
        cusplat_sh_amie = 75.0
        cuspmlt_sh_amie = 11.
      endif
      if (cusplat_sh_amie < 60.0) then
        cusplat_sh_amie = 60.0
        cuspmlt_sh_amie = 11.
      endif
!     cusplat_sh_amie = amin1(65.0,cusplat_sh_amie)

!     if (cusplat_nh_amie > 65.0) then
!       cusplat_nh_amie = 65.0
!       cuspmlt_nh_amie = 11.
!     endif
      if (cusplat_nh_amie > 75.0) then
        cusplat_nh_amie = 75.0
        cuspmlt_nh_amie = 11.
      endif
      if (cusplat_nh_amie < 60.0) then
        cusplat_nh_amie = 60.0
        cuspmlt_nh_amie = 11.
      endif
!     cusplat_nh_amie = amin1(65.0,cusplat_nh_amie)
      if (cuspmlt_sh_amie > 12.5) cuspmlt_sh_amie = 12.5
      if (cuspmlt_sh_amie < 11.0) cuspmlt_sh_amie = 11.0
      if (cuspmlt_nh_amie > 12.5) cuspmlt_nh_amie = 12.5
      if (cuspmlt_nh_amie < 11.0) cuspmlt_nh_amie = 11.0
      crad(1) = (90.-cusplat_sh_amie)*pi/180.
      crad(2) = (90.-cusplat_nh_amie)*pi/180.
      phida(1) = (cuspmlt_sh_amie - 12.) * pi / 12.
      phida(2) = (cuspmlt_nh_amie - 12.) * pi / 12.
!
! gl - 14/07/2002 changed from Emery's original amieterp.f
!
!  Change mag coords so latitude extends to the equator in into the opposite
!   pole, and so MLT is converted to apex longitues between -180 and +180.
!  Put EPOTDUP, EMKEV and WPM2 into S. Hem and N. Hem alike in EPOTM, EKEVM,
!   and EFLXM, dimensioned (IMXMP,JMNH) where MLAT goes from -90 to ~-2
!   for ihsoln=1, and from +90 to ~+2 for ihsoln=2 since GRDSET is for
!   -90 to the equator.  (Will work for +90 to the equator also.)
!   AMIE grid goes from -90 to ~-40 for SH and +90 to ~+40 for NH.
!   TGCM grid goes from -87.5 to -2.5 for SH and +2.5 to +87.5 for NH.
!   Hence, at very end, need to reverse order of NH in whole globe arrays.
!
!     do i=1,lonp1
!  Initialize arrays around equator
!       do j=latp1+1,ithmx
!         potm(i,j) = 0.
!         potm(i,jmxm+1-j) = 0.
!         ekvm(i,j) = ekv_sh_amie(i,latp1)
!         ekvm(i,jmxm+1-j) = ekv_nh_amie(i,latp1)
!         efxm(i,j) = 0.
!         efxm(i,jmxm+1-j) = 0.
!       enddo
!  Put in AMIE arrays from pole to JMFT
!       do j=1,latp1
!         potm(i,j) = pot_sh_amie(i,j)
!         potm(i,latp1+1-j) = pot_nh_amie(i,j)
!         ekvm(i,j) = ekv_sh_amie(i,j)
!         ekvm(i,latp1+1-j) = ekv_nh_amie(i,j)
!         efxm(i,j) = efx_sh_amie(i,j)
!         efxm(i,latp1+1-j) = efx_nh_amie(i,j)
!       enddo
!     enddo

! mlongitude starts from 180 degree
      rot = 250./15. - iutsec/3600.
      dmltm = 24./float(lonmx)
      do i=1,lonp1
        xmlt = (i-1)*dmltm - rot + 24.
        xmlt = AMOD(xmlt,24.)
        m = IFIX(xmlt/dmltm + 1.01)
        mp1 = m + 1
        IF (mp1 > lonp1) mp1 = 2
        del = xmlt - (m-1)*dmltm
C  Initialize arrays around equator
        do j=latp1+1,ithmx
          potm(i,j) = 0.
          potm(i,jmxm+1-j) = 0.
          ekvm(i,j) = (1.-del)*ekv_sh_amie(m,latp1) + 
     |                        del*ekv_sh_amie(mp1,latp1)
          ekvm(i,jmxm+1-j) = (1.-del)*ekv_nh_amie(m,latp1) +
     |                        del*ekv_nh_amie(mp1,latp1)
          efxm(i,j) = 0.
          efxm(i,jmxm+1-j) = 0.
        enddo
C  Put in AMIE arrays from pole to latp1
        do j=1,latp1
          potm(i,j) = (1.-del)*pot_sh_amie(m,j) + del*pot_sh_amie(mp1,j)
          potm(i,jmxm+1-j) = (1.-del)*pot_nh_amie(m,j) + 
     |                       del*pot_nh_amie(mp1,j)
          ekvm(i,j) = (1.-del)*ekv_sh_amie(m,j) + del*ekv_sh_amie(mp1,j)
          ekvm(i,jmxm+1-j) = (1.-del)*ekv_nh_amie(m,j) + 
     |                       del*ekv_nh_amie(mp1,j)
          efxm(i,j) = (1.-del)*efx_sh_amie(m,j) + del*efx_sh_amie(mp1,j)
          efxm(i,jmxm+1-j) = (1.-del)*efx_nh_amie(m,j) + 
     |                       del*efx_nh_amie(mp1,j)
        enddo
     
      enddo

C  Set up coeffs to go between EPOTM(IMXMP,JMNH) and TIEPOT(IMAXM,JMAXMH)
C
C     ****     SET GRID SPACING DLATM, DLONG, DLONM
C  DMLAT=lat spacing in degrees of AMIE apex grid
      dtr = pi/180.
      DMLAT = 180. / FLOAT(jmxm-1)
      DLATM = DMLAT*dtr
      DLONM = 2.*PI/FLOAT(lonmx)
      DMLTM = 24./FLOAT(lonmx)
C     ****
C     ****     SET ARRAY YLATM (LATITUDE VALUES FOR GEOMAGNETIC GRID
C     ****
      alatm(1) = -PI/2.
      alat(1) = -90.
      alatm(jmxm) = PI/2.
      alat(jmxm) = 90.
      do i = 2,ithmx
        alat(i) = alat(i-1)+dlatm/dtr
        alat(jmxm+1-i) = alat(jmxm+2-i)-dlatm/dtr
        alatm(i) = alatm(i-1)+dlatm
        alatm(jmxm+1-i) = alatm(jmxm+2-i)-dlatm
      enddo
      alon(1) = -pi/dtr
      alonm(1) = -pi
      do i=2,lonp1
        alon(i) = alon(i-1) + dlonm/dtr
        alonm(i) = alonm(i-1) + dlonm
      enddo

! ylatm and ylonm are arrays of latitudes and longitudes of the
! distored magnetic grids in radian - from consdyn.h
!  Convert from apex magnetic grid to distorted magnetic grid
!
! Allocate workspace for regrid routine rgrd2.f:
      lw = nmlonp1+nmlat+2*nmlonp1
      if (.not. allocated(w)) allocate(w(lw),stat=ier)
      if (ier /= 0) write(6,"('>>> horizontal_interp: error allocating',
     |  ' w(lw): lw=',i6,' ier=',i4)") lw,ier
      liw = nmlonp1 + nmlat
      if (.not. allocated(iw)) allocate(iw(liw),stat=ier)
      if (ier /= 0) write(6,"('>>> horzontal_interp: error allocating',
     |  ' iw(liw): liw=',i6,' ier=',i4)") liw,ier
      intpol(:) = 1 ! linear (not cubic) interp in both dimensions
      if (alatm(1) > ylatm(1)) alatm(1) = ylatm(1)
      if (alatm(jmxm) < ylatm(nmlat)) alatm(jmxm) = ylatm(nmlat)
      if (alonm(1) > ylonm(1)) alonm(1) = ylonm(1)
      if (alonm(lonp1) < ylonm(nmlonp1)) alonm(lonp1) = ylonm(nmlonp1)
      call rgrd2(lonp1,jmxm,alonm,alatm,potm,nmlonp1,nmlat,
     |           ylonm,ylatm,tiepot,intpol,w,lw,iw,liw,ier)
      call rgrd2(lonp1,jmxm,alonm,alatm,ekvm,nmlonp1,nmlat,
     |           ylonm,ylatm,tieekv,intpol,w,lw,iw,liw,ier)
      call rgrd2(lonp1,jmxm,alonm,alatm,efxm,nmlonp1,nmlat,
     |           ylonm,ylatm,tieefx,intpol,w,lw,iw,liw,ier)

C  Convert from TGCM distorted magnetic grid to geographic one
!      CALL GRDTRP(tiepot(1,1),potg(1,0),im(1,0),jm(1,0),
!    |    dim(1,0),djm(1,0),IMAXGP,IMAXMP,IMAXG,JMAXG+2,IMAXM,JMAXM)
!      CALL GRDTRP(tieekv(1,1),ekvg(1,0),im(1,0),jm(1,0),
!    |    DIM(1,0),DJM(1,0),IMAXGP,IMAXMP,IMAXG,JMAXG+2,IMAXM,JMAXM)
!      CALL GRDTRP(tieefx(1,1),efxg(1,0),im(1,0),jm(1,0),
!    |    DIM(1,0),DJM(1,0),IMAXGP,IMAXMP,IMAXG,JMAXG+2,IMAXM,JMAXM)

       call mag2geo(tiepot(1,1),potg(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
       call mag2geo(tieekv(1,1),ekvg(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
       call mag2geo(tieefx(1,1),efxg(1,0),im(1,0),jm(1,0),
     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
C     ****
C     ****     INSERT PERIODIC POINTS
C     ****
      DO j = 1,nlat
        ekvg(nlonp1,j) = ekvg(1,j)
        efxg(nlonp1,j) = efxg(1,j)
        potg(nlonp1,j) = potg(1,j)
      ENDDO

! temp - save the AMIE output in geographic grid
      amie_debug = 1
      if (amie_debug > 0) then
!     do j=1,nmlat
!       do i=1,nmlonp1
!         tiepot_sech(i,:) = tiepot(i,j)
!       enddo
!       call addfld('TIEPOT',' ',' ',tiepot_sech(:,1:nlevp1),
!    |    'mlon',1,nmlonp1,'mlev',1,nlevp1,j)
!     enddo
! tiepot(nmlonp1,nmlat)
!       call addfsech_ij('TIEPOT',' ',' ',tiepot,1,nmlonp1,1,nmlat)
!       call addfsech_ij('TIEEKV',' ',' ',tieekv,1,nmlonp1,1,nmlat)
!       call addfsech_ij('TIEEFX',' ',' ',tieefx,1,nmlonp1,1,nmlat)

        call addfld('TIEPOT',' ',' ',tiepot,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
        call addfld('TIEEKV',' ',' ',tieekv,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)
        call addfld('TIEEFX',' ',' ',tieefx,'mlon',1,nmlonp1,
     |    'mlat',1,nmlat,0)

        do i=1,nlon
          potg_sech(i+2,:) = potg(i,1:nlat)
          ekvg_sech(i+2,:) = ekvg(i,1:nlat)
          efxg_sech(i+2,:) = efxg(i,1:nlat)
        enddo
        do i=1,2
          potg_sech(i,:) = potg_sech(nlon+i,:)
          potg_sech(nlon+i+2,:) = potg_sech(i,:)
          ekvg_sech(i,:) = ekvg_sech(nlon+i,:)
          ekvg_sech(nlon+i+2,:) = ekvg_sech(i,:)
          efxg_sech(i,:) = efxg_sech(nlon+i,:)
          efxg_sech(nlon+i+2,:) = efxg_sech(i,:)
        enddo
!       call addfsech_ij('POTG',' ',' ',potg_sech,1,nlonp4,1,nlat)
!       call addfsech_ij('EKVG',' ',' ',ekvg_sech,1,nlonp4,1,nlat)
!       call addfsech_ij('EFXG',' ',' ',efxg_sech,1,nlonp4,1,nlat)

        call addfld('POTG',' ',' ',potg_sech,'lon',1,nlonp4,
     |    'lat',1,nlat,0)
        call addfld('EKVG',' ',' ',ekvg_sech,'lon',1,nlonp4,
     |    'lat',1,nlat,0)
        call addfld('EFXG',' ',' ',efxg_sech,'lon',1,nlonp4,
     |    'lat',1,nlat,0)
      endif
!      
      if (iprint > 0) then
        write(6,"('getamie: AMIE data interpolated to date and time')")
        write(6,"(72('-'),/)")
      endif
      end subroutine getamie    
!-------------------------------------------------------------------
      subroutine boxcar_ave(x,y,lon,lat,mtime,itime,ibox)
!
! perform boxcar average
!
! Args:
      integer,intent(in) :: lon,lat,mtime,itime,ibox
      real,dimension(lon,lat,mtime), intent(in) :: x
      real,dimension(lon,lat), intent(out) :: y
! Local:
      integer :: i, iset, iset1
!
      iset = itime - ibox/2
      if (iset < 1) iset = 1
      iset1 = iset + ibox
      if (iset1 > mtime) then
        iset1 = mtime
        iset = iset1 - ibox
      endif
!     write(6,"('boxcar_ave: mtime,itime,ibox',3i5)")
!    |  mtime,itime,ibox
!
      y(:,:) = 0.
      do i=iset,iset1
        y(:,:) = y(:,:) + x(:,:,i)
      enddo     
      if (ibox > 0) y(:,:) = y(:,:)/ibox
!
      end subroutine boxcar_ave
!-----------------------------------------------------------------------
      subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg,nlonm,
     |  nlatm)
!
! Args:
      integer,intent(in) :: lg,lm,nlong,nlatg,nlonm,nlatm
      integer,intent(in) :: im(lg,*),jm(lg,*)
      real,intent(in) :: am(lm,*),dim(lg,*),djm(lg,*)
      real,intent(out) :: ag(lg,*)
!
! Local:
      integer :: ig,jg
!
      do jg=1,nlatg 
        do ig=1,nlong
          ag(ig,jg) =
     |    am(im(ig,jg)  ,jm(ig,jg))  *(1.-dim(ig,jg))*(1.-djm(ig,jg))+
     |    am(im(ig,jg)+1,jm(ig,jg))  *    dim(ig,jg) *(1.-djm(ig,jg))+
     |    am(im(ig,jg)  ,jm(ig,jg)+1)*(1.-dim(ig,jg))*djm(ig,jg)+
     |    am(im(ig,jg)+1,jm(ig,jg)+1)*    dim(ig,jg) *djm(ig,jg)
        enddo ! ig=1,nlong
      enddo ! jg=1,nlatg 
      end subroutine mag2geo
!-----------------------------------------------------------------------
      subroutine rpt_ncerr(istat,msg)
!
! Handle a netcdf lib error:
!
      integer,intent(in) :: istat
      character(len=*),intent(in) :: msg
!
      write(6,"(/72('-'))")
      write(6,"('>>> Error from netcdf library:')")
      write(6,"(a)") trim(msg)
      write(6,"('istat=',i5)") istat
      write(6,"(a)") nf_strerror(istat)
      write(6,"(72('-')/)")
      return
      end subroutine rpt_ncerr
      end module amie_module
