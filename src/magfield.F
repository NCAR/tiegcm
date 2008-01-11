!
      module magfield_module
!
! Read and calculate magnetic field data used for transformation 
!   between geographic and geomagnetic coordinate systems. 
!   Sub magfield is called once per run from main tgcm.F.
!   Magfield calls nc_rdmag, magdyn, and dynpotmag.
!
! Nc_rdmag:  Read netcdf file containing magnetic field information.
! Magdyn:    Establish magnetic field and grid (bx,by,bz).
! Dynpotmag: Transform electric potential (from source history) from
!            geographic to geomagnetic grid.
!
! Sub sunloc is called once per timestep from advance to determine
!   sun's longitudes for current ut model time.
!
      use params_module,only: nmlonp1,nmlat,nlat,nlonp4,nlatp2,nlonp2,
     |  nlon,nlonp1,nlevp1,nmlon,nlatp1,nmlev
      use input_module,only: tempdir,magvol
      use nchist_module,only: nc_open,nc_close,handle_ncerr
      implicit none
!
! Netcdf header file:
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
!
! Grid dimensions:
!
! Magnetic quantities read from magvol by sub nc_rdmag (called by rdmag)
! (see input_module for magvol file) (these were in old version fieldz.h)
!
      real,dimension(nlonp1,0:nlatp1) ::
     |  alatm,   ! geomagnetic latitude at each geographic grid point (radians)
     |  alonm,   ! geomagnetic longitude at each geographic grid point (radians)
     |  xb,      ! northward component of magnetic field
     |  yb,      ! eastward component of magnetic field
     |  zb,      ! downward component of magnetic field
     |  bmod,    ! magnitude of magnetic field (gauss?)
     |  dmlat,   ! dipole latitude corresponding to apex of field line
     |  rmag11,  ! (a1.a1)/p*sin(i)*cos(thetas)
     |  rmagc,   ! (a1.a2)/p*sin(i)
     |  rmag2,   ! 1./bmod
     |  rmag22,  ! (a2.a2)/p*sin(i)/cos(thetas)
     |  rjacd,   ! determinant of rjac
     |  p
!
! rjac: scaled derivatives of geomagnetic coords wrt geographic coordinates.
! rjac(1,1) = cos(thetas)/cos(theta)*d(lamdas)/d(lamda)
! rjac(1,2) = cos(thetas)*d(lamdas)/d(theta)
! rjac(2,1) = 1./cos(theta)*d(thetas)/d(lamda)
! rjac(2,2) = d(thetas)/d(theta)
! where (lamda,theta) are geographic coordinates
!       (lamdas,thetas) are geomagnetic coordinates
!
      real :: rjac(nlonp1,0:nlatp1,2,2) 
!
! av  = the two magnetic vectors a1 and a2
! av1 = a1
! av2 = a2/cos(thetas)
!
      real :: av  (nlonp1,0:nlatp1,3,2)
!
! Quantities needed to transform scalar fields between geographic and
! geomagnetic coordinate systems (these were in old version cterp.h)
!
      integer :: 
     |  ig(nmlonp1,nmlat),   ! geog lon grid containing each geomag point
     |  jg(nmlonp1,nmlat),   ! geog lat grid containing each geomag point
     |  im(nlonp1,0:nlatp1), ! geomag lon grid containing each geog point
     |  jm(nlonp1,0:nlatp1)  ! geomag lat grid containing each geog point
!
! wt(4) are interpolation weights to be applied to function values at 4 
!   corners of geographic grid element (ig,jg)
! dim and djm are fractions in i and j directions, used for bilinear 
!   interpolation in geomagnetic grid element.
!
      real ::
     |  wt(4,nmlonp1,nmlat),   
     |  dim(nlonp1,0:nlatp1),
     |  djm(nlonp1,0:nlatp1)
!
! Trigonometric factors needed in the calculation of the derivatives
! of the geomagnetic coordinates wrt the geographic coordinates
! (these were in old versions header file trig.h) 
!
      real,dimension(nlonp1,nlat) ::
     |  cslatm, ! cos(thetas)
     |  snlatm, ! sin(thetas)
     |  cslonm, ! cos(lamdas)
     |  snlonm  ! sin(lamdas)
!
! 10/17/02 bf: cslatg and snlatg dim changed from nlat to 0:nlatp1
!   for apxparm (apex.F).
      real ::
     |  cslatg(0:nlatp1), ! cos(theta)
     |  snlatg(0:nlatp1), ! sin(theta)
     |  cslong(nlonp1),   ! cos(lamda)
     |  snlong(nlonp1)    ! sin(lamda)
!
! Magnetic field at geographic grid. These are defined by magdyn after 
!   reading magdat file (formerly in common /MAGFLD/ and header file trgm.h)
! (note bx,by,bz,bmod2 are dimensioned nlonp4 rather than -1:nlonp2, as
!  in earlier versions, so these fields can be referenced at i rather than
!  i-2, e.g. in sub oplus)
!
      real,dimension(nlonp4,-1:nlatp2) :: bx,by,bz,bmod2
      real,dimension(nlonp4,nlat) ::
     |  rlatm,rlonm,dipmag,decmag,sndec,csdec,sn2dec,sncsdc
!
! sunlons: sun's longitude in dipole coordinates (see sub sunloc)
! (this was dlons in earlier versions)
!
      real :: sunlons(nlat)
!
      contains
!-----------------------------------------------------------------------
      subroutine magfield(dynamo)
!
! Read magnetic field data file:
!
      implicit none
!
! Args:
      integer,intent(in) :: dynamo
!
! Local:
      integer :: nwds,j,k,i,ii,ier
      character(len=80) :: dskfile
!
! Get and read magnetic data file.
! 2/00: using new netcdf mag file (see ~foster/tgcm/mkmag)
! 9/00: In new getms, if MPI, only master proc actually acquires mspath.
! For systems without ncar mss (MSS==0), MAGVOL given in user input file 
!   should be a local file name, e.g., MAGVOL='TGCM.data.magfield.nc'

! btf 1/14/05: Removing calls to acquire and read magvol, since apex coord
!              code is always called if DYNAMO=1. (i.e., there is no
!              "old dynamo", and the "new dynamo" always uses apex code.
!
#if (MSS == 0)
      dskfile = magvol
      call get_diskfile(dskfile,' ') 
#else
      dskfile = ' '
      call getfile(magvol,dskfile)
#endif
!
! If dynamo==1, then dynamo will be called. The dynamo uses
! apxparm (apex.F), so nc_rdmag is not necessary in this case.
!
      if (dynamo <= 0) then
        call nc_rdmag(dskfile)
      else
        write(6,"('magfield: nc_rdmag not called because dynamo=',
     |    i2)") dynamo
      endif
!
      call magdyn
      if (dynamo > 0) call dynpotmag
      end subroutine magfield
!-------------------------------------------------------------------
      subroutine nc_rdmag(dskfile)
! 
! Read netcdf magnetic data file.
! The mss data file (magvol) is set from input according to current
!   horizontal resolution. Currently (1/02) there are 2 such files:
!   /TGCM/data/magdat.nc for 5.0 degree resolution, and
!   /TGCM/data/magdat_2.5h.nc for 2.5 degree resolution.
! These files are generated by code in hao /home/tgcm/mkmag, using
!   the apex parm code (s.a. Roy Barnes).
!
! The "original" cray-blocked file for 5 degrees was
!   /ECRIDLEY/ECR90/ECRMG6 (s.a. code in ~foster/tgcm/mkmag)
! 
      implicit none
! Args:
      character(len=*),intent(in) :: dskfile
!
! Local:
      integer :: istat,ncid,j
      integer :: ids1(1),ids2(2),ids3(3),ids4(4)
      integer :: id_nlat,id_nlonp1,id_nlatp2,id_dim2,id_dim3,
     |  id_nmlonp1,id_nmlat,id_dim4
      integer :: idv_alatm,idv_alonm,idv_xb,idv_yb,idv_zb,idv_bmod,
     |  idv_dmlat,idv_rjac,idv_av,idv_p,idv_rmag11,idv_rmagc,idv_rmag2,
     |  idv_rmag22,idv_rjacd,idv_im,idv_jm
      integer :: idv_ig,idv_jg,idv_wt,idv_dim,idv_djm
      integer :: idv_cslatm,idv_snlatm,idv_cslonm,idv_snlonm,idv_cslatg,
     |  idv_snlatg,idv_cslong,idv_snlong
      integer :: start_1d(1),count_1d(1),start_2d(2),count_2d(2),
     |  start_3d(3),count_3d(3),start_4d(4),count_4d(4)
      character(len=8) :: dimname
      real :: fmin,fmax
!
      write(6,"(/,72('-'))") 
      write(6,"('nc_rdmag: Read magnetic field data file ',a,
     |  /,10x,' (mss file ',a,')')") trim(dskfile),trim(magvol)
!
! Open the netcdf dataset:
      call nc_open(ncid,dskfile,'OLD','READ')
      if (ncid <= 0) then
        write(6,"(/,'>>> nc_rdmag: error opening netcdf mag data ',
     |    'file ',a)") trim(dskfile)
        call shutdown('nc_rdmag')
!     else
!       write(6,"('nc_rdmag: opened netcdf mag data file ',a,
!    |    ' ncid=',i8)") trim(dskfile),ncid
      endif 
!
! Check dimensions:
      call checkdim(ncid,"nlat"   ,nlat)
      call checkdim(ncid,"nlonp1" ,nlonp1)
      call checkdim(ncid,"nlatp2" ,nlatp1+1)
!     call checkdim(ncid,"nmlonp1",nmlonp1)
      call checkdim(ncid,"nmlon"  ,nmlonp1)
      call checkdim(ncid,"nmlat"  ,nmlat)
!
! Read variables for fieldz.h:
! 2-d doubles (nlonp1,0:nlatp1):
      start_2d(:) = 1
      count_2d(1) = nlonp1
      count_2d(2) = nlatp2
      call rd2dfld(ncid,'ALATM ',idv_alatm ,start_2d,count_2d,alatm )
      call rd2dfld(ncid,'ALONM ',idv_alonm ,start_2d,count_2d,alonm )
      call rd2dfld(ncid,'XB    ',idv_xb    ,start_2d,count_2d,xb    )
      call rd2dfld(ncid,'YB    ',idv_yb    ,start_2d,count_2d,yb    )
      call rd2dfld(ncid,'ZB    ',idv_zb    ,start_2d,count_2d,zb    )
      call rd2dfld(ncid,'BMOD  ',idv_bmod  ,start_2d,count_2d,bmod  )
      call rd2dfld(ncid,'DMLAT ',idv_dmlat ,start_2d,count_2d,dmlat )
      call rd2dfld(ncid,'P     ',idv_p     ,start_2d,count_2d,p     )
      call rd2dfld(ncid,'RMAG11',idv_rmag11,start_2d,count_2d,rmag11)
      call rd2dfld(ncid,'RMAGC ',idv_rmagc ,start_2d,count_2d,rmagc )
      call rd2dfld(ncid,'RMAG2 ',idv_rmag2 ,start_2d,count_2d,rmag2 )
      call rd2dfld(ncid,'RMAG22',idv_rmag22,start_2d,count_2d,rmag22)
      call rd2dfld(ncid,'RJACD ',idv_rjacd ,start_2d,count_2d,rjacd )
!
! RJAC(nlonp1,0:nlatp1,2,2):
      start_4d(:) = 1
      count_4d(1) = nlonp1
      count_4d(2) = nlatp2
      count_4d(3:4) = 2
      istat = nf_inq_varid(ncid,'RJAC',idv_rjac)
      istat = nf_get_vara_double(ncid,idv_rjac,start_4d,count_4d,rjac)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for rjac')
      call fminmax(rjac,nlonp1*nlatp2*2*2,fmin,fmax)
      write(6,"('  RJAC  min,max=',2e12.4)") fmin,fmax
!
! AV(nlonp1,0:nlatp1,3,2):
      start_4d(:) = 1
      count_4d(1) = nlonp1
      count_4d(2) = nlatp2
      count_4d(3) = 3
      count_4d(4) = 2
      istat = nf_inq_varid(ncid,'AV',idv_av)
      istat = nf_get_vara_double(ncid,idv_av,start_4d,count_4d,av)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for av')
      call fminmax(av,nlonp1*nlatp2*3*2,fmin,fmax)
      write(6,"('  AV    min,max=',2e12.4)") fmin,fmax
!
! Read variables to module data:
!
! IG(nmlonp1,nmlat)
      start_2d(:) = 1
      count_2d(1) = nmlonp1
      count_2d(2) = nmlat
      istat = nf_inq_varid(ncid,'IG',idv_ig)
      istat = nf_get_vara_int(ncid,idv_ig,start_2d,count_2d,ig)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_int for ig')
!
! JG(nmlonp1,nmlat)
      count_2d(1) = nmlonp1
      count_2d(2) = nmlat
      istat = nf_inq_varid(ncid,'JG',idv_jg)
      istat = nf_get_vara_int(ncid,idv_jg,start_2d,count_2d,jg)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_int for jg')
!
! IM(nlonp1,0:nlatp1):
      count_2d(1) = nlonp1
      count_2d(2) = nlatp2
      istat = nf_inq_varid(ncid,'IM',idv_im)
      istat = nf_get_vara_int(ncid,idv_im,start_2d,count_2d,im)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_int for im')
!
! JM(nlonp1,0:nlatp1):
      istat = nf_inq_varid(ncid,'JM',idv_jm)
      istat = nf_get_vara_int(ncid,idv_jm,start_2d,count_2d,jm)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_int for jm')
!
! WT(4,nmlonp1,nmlat):
      start_3d(:) = 1
      count_3d(1) = 4
      count_3d(2) = nmlonp1
      count_3d(3) = nmlat
      istat = nf_inq_varid(ncid,'WT',idv_wt)
      istat = nf_get_vara_double(ncid,idv_wt,start_3d,count_3d,wt)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for wt')
      call fminmax(wt,4*nmlonp1*nmlat,fmin,fmax)
      write(6,"('  WT    min,max=',2e12.4)") fmin,fmax
!
! DIM and DJM(nlonp1,0:nlatp1):
      start_2d(:) = 1
      count_2d(1) = nlonp1
      count_2d(2) = nlatp2
      call rd2dfld(ncid,'DIM   ',idv_dim ,start_2d,count_2d,dim )
      call rd2dfld(ncid,'DJM   ',idv_djm ,start_2d,count_2d,djm )
!
! cslatm, snlatm, cslonm, snlonm:
      start_2d(:) = 1
      count_2d(1) = nlonp1
      count_2d(2) = nlat
      call rd2dfld(ncid,'CSLATM ',idv_cslatm ,start_2d,count_2d,cslatm)
      call rd2dfld(ncid,'SNLATM ',idv_snlatm ,start_2d,count_2d,snlatm)
      call rd2dfld(ncid,'CSLONM ',idv_cslonm ,start_2d,count_2d,cslonm)
      call rd2dfld(ncid,'SNLONM ',idv_snlonm ,start_2d,count_2d,snlonm)
!
! CSLATG(nlat):
      start_1d(:) = 1
      count_1d(1) = nlat
      istat = nf_inq_varid(ncid,'CSLATG',idv_cslatg)
      istat = nf_get_vara_double(ncid,idv_cslatg,start_1d,count_1d,
     |  cslatg)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for cslatg')
      call fminmax(cslatg,nlat,fmin,fmax)
      write(6,"('  CSLATG  min,max=',2e12.4)") fmin,fmax
!
! SNLATG(nlat):
      istat = nf_inq_varid(ncid,'SNLATG',idv_snlatg)
      istat = nf_get_vara_double(ncid,idv_snlatg,start_1d,count_1d,
     |  snlatg)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for snlatg')
      call fminmax(snlatg,nlat,fmin,fmax)
      write(6,"('  SNLATG  min,max=',2e12.4)") fmin,fmax
!
! CSLONG(nlonp1):
      count_1d(1) = nlonp1
      istat = nf_inq_varid(ncid,'CSLONG',idv_cslong)
      istat = nf_get_vara_double(ncid,idv_cslong,start_1d,count_1d,
     |  cslong)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for cslong')
      call fminmax(cslong,nlonp1,fmin,fmax)
      write(6,"('  CSLONG  min,max=',2e12.4)") fmin,fmax
!
! SNLONG(nlonp1):
      istat = nf_inq_varid(ncid,'SNLONG',idv_snlong)
      istat = nf_get_vara_double(ncid,idv_cslong,start_1d,count_1d,
     |  snlong)
      if (istat /= NF_NOERR) call handle_ncerr(istat,
     |    'Error return from nf_get_vara_double for snlong')
      call fminmax(snlong,nlonp1,fmin,fmax)
      write(6,"('  SNLONG  min,max=',2e12.4)") fmin,fmax
!
! Close the dataset:
      call nc_close(ncid)
      write(6,"('Completed read of magnetic field data file.')")
      write(6,"(72('-'),/)")
      end subroutine nc_rdmag
!-------------------------------------------------------------------
      subroutine checkdim(ncid,dimname,iparam)
!
! Get length of dimension "dimname". If this length is not equal
! to iparam, stop with error message.
!
      implicit none
!
! Args:
      integer,intent(in) :: ncid,iparam
      character(len=*),intent(in) :: dimname
!
! Local:
      integer :: istat,iddim,len
      character(len=80) :: char80
!
! Get dim id:
      istat = nf_inq_dimid(ncid,dimname,iddim)
      if (istat /= NF_NOERR) then
        write(char80,"('nc_rdmag: error getting dim id for ',
     |    a)") dimname
        call handle_ncerr(istat,char80)
      endif
!
! Get dim length:
      istat = nf_inq_dimlen(ncid,iddim,len)
      if (istat /= NF_NOERR) then
        write(char80,"('nc_rdmag: error getting length of ',
     |    'dimension ',a)") dimname
        call handle_ncerr(istat,char80)
      endif
!
! Compare with iparam:
      if (len /= iparam) then
        write(6,"(/,'>>> nc_rdmag: unexpected length for ',
     |    'dimension ',a)") dimname
        write(6,"('  length read = ',i3,' should be = ',i3)") 
     |    len,iparam
        call shutdown('nc_rdmag')
      endif
      end subroutine checkdim
!-------------------------------------------------------------------
      subroutine rd2dfld(ncid,name,idvout,start_2d,count_2d,var)
      implicit none
!
! Read 2-d double array from ncid to var:
!
! Args:
      integer,intent(in) :: ncid,start_2d(2),count_2d(2)
      character(len=*),intent(in) :: name
      integer,intent(out) :: idvout
      real,intent(out) :: var(count_2d(1),count_2d(2))
!
! Local:
      integer :: istat
      character(len=80) :: char80
      real :: fmin,fmax
!
      istat = nf_inq_varid(ncid,name,idvout)
      istat = nf_get_vara_double(ncid,idvout,start_2d,count_2d,var)
      write(char80,"('Error return from nf_get_vara_double for var',
     |  a)") name
      if (istat /= NF_NOERR) call handle_ncerr(istat,char80)
!     call fminmax(var,count_2d(1)*count_2d(2),fmin,fmax)
!     write(6,"('rd2dfld: ',a,' fmin,max=',2e12.4)") name,fmin,fmax
      end subroutine rd2dfld
!-----------------------------------------------------------------------
#include <defs.h>
      subroutine magdyn
!
! Local:
!     real,parameter :: sin10=0.17   ! should be equal to dipmin
      real :: cos10
!
! sin10 should be the same as dipmin (see cons.F)
!
#if (NLAT==36 && NLON==72)
      real,parameter :: sin10=0.17   ! 5.0 deg horizontal resolution
#elif (NLAT==72 && NLON==144)
      real,parameter :: sin10=0.24   ! 2.5 deg horizontal resolution
#else
      UNKNOWN NLAT,NLON      ! compilation will stop here if unknown res
#endif
!
! Field are calculated in local arrays dimensioned -1:nlonp2, then transferred 
!   to module data bx,by,bz,bmod2 at 1:nlonp4 for use in the model.
!
      real,dimension(-1:nlonp2,-1:nlatp2) :: bxtmp,bytmp,bztmp,bmod2tmp
      integer :: i,j
!
      cos10 = sqrt(1.-sin10**2)
      do j = 1,nlat
        do i = 1,nlon
          rlatm(i+2,j) = alatm(i,j)
          rlonm(i+2,j) = alonm(i,j)
          dipmag(i+2,j) = atan(zb(i,j)/sqrt(xb(i,j)**2+yb(i,j)**2))
          decmag(i+2,j) = atan2(yb(i,j),xb(i,j))
          sndec(i+2,j) = sin(decmag(i+2,j))
          csdec(i+2,j) = cos(decmag(i+2,j))
          sn2dec(i+2,j) = sndec(i+2,j)**2
          sncsdc(i+2,j) = sndec(i+2,j)*csdec(i+2,j)
          bxtmp(i,j) = yb(i,j)/bmod(i,j)
          bytmp(i,j) = xb(i,j)/bmod(i,j)
          bztmp(i,j) = -zb(i,j)/bmod(i,j)
          bmod2tmp(i,j) = bmod(i,j)
!
! Set minimum dip to 10 degrees
          if (abs(bztmp(i,j))-sin10 < 0.) then
            bxtmp(i,j)=bxtmp(i,j)*(cos10/sqrt(1.-bztmp(i,j)**2))
            bytmp(i,j)=bytmp(i,j)*(cos10/sqrt(1.-bztmp(i,j)**2))
            bztmp(i,j)=sign(sin10,bztmp(i,j))
          endif
        enddo
      enddo
!
! Values at j = -1, 0, nlatp1, nlatp2:
      do j = 1,2
        do i = 1,nlon
          bxtmp(i,j-2) = -bxtmp(1+mod(i-1+nlon/2,nlon),3-j)
          bytmp(i,j-2) = -bytmp(1+mod(i-1+nlon/2,nlon),3-j)
          bztmp(i,j-2) =  bztmp(1+mod(i-1+nlon/2,nlon),3-j)
          bmod2tmp(i,j-2) = bmod2tmp(1+mod(i-1+nlon/2,nlon),3-j)
          bxtmp(i,nlat+j) = -bxtmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bytmp(i,nlat+j) = -bytmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bztmp(i,nlat+j) =  bztmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
          bmod2tmp(i,nlat+j) = bmod2tmp(1+mod(i-1+nlon/2,nlon),nlat+1-j)
        enddo
      enddo
!
! Periodic points:
      do i = 1,2
        do j = 1,nlat
          rlatm (i,j) = rlatm (i+nlon,j)
          rlonm (i,j) = rlonm (i+nlon,j)
          dipmag(i,j) = dipmag(i+nlon,j)
          decmag(i,j) = decmag(i+nlon,j)
          sndec (i,j) = sndec (i+nlon,j)
          csdec (i,j) = csdec (i+nlon,j)
          sn2dec(i,j) = sn2dec(i+nlon,j)
          sncsdc(i,j) = sncsdc(i+nlon,j)
!
          rlatm (i+nlonp2,j) = rlatm (i+2,j)
          rlonm (i+nlonp2,j) = rlonm (i+2,j)
          dipmag(i+nlonp2,j) = dipmag(i+2,j)
          decmag(i+nlonp2,j) = decmag(i+2,j)
          sndec (i+nlonp2,j) = sndec (i+2,j)
          csdec (i+nlonp2,j) = csdec (i+2,j)
          sn2dec(i+nlonp2,j) = sn2dec(i+2,j)
          sncsdc(i+nlonp2,j) = sncsdc(i+2,j)
        enddo
      enddo
      do i = 1,2
        do j = -1,nlatp2
          bxtmp(i-2,j) = bxtmp(i-2+nlon,j)
          bytmp(i-2,j) = bytmp(i-2+nlon,j)
          bztmp(i-2,j) = bztmp(i-2+nlon,j)
          bmod2tmp(i-2,j) = bmod2tmp(i-2+nlon,j)

          bxtmp(i+nlon,j) = bxtmp(i,j)
          bytmp(i+nlon,j) = bytmp(i,j)
          bztmp(i+nlon,j) = bztmp(i,j)
          bmod2tmp(i+nlon,j) = bmod2tmp(i,j)
        enddo
      enddo
!
! Save local fields -1:nlonp2 to module data at 1:nlonp4.
! This way, model routines (e.g., oplus) can reference i rather than i-2.
      do i=1,nlonp4
        bx(i,:) = bxtmp(i-2,:)
        by(i,:) = bytmp(i-2,:)
        bz(i,:) = bztmp(i-2,:)
        bmod2(i,:) = bmod2tmp(i-2,:)
      enddo ! i=1,nlonp4

!     do j=-1,nlatp2
!       write(6,"('magdyn: j=',i2)") j
!       write(6,"('bx(:,j)=',/,(6e12.4))") bx(:,j)
!       write(6,"('by(:,j)=',/,(6e12.4))") by(:,j)
!       write(6,"('bz(:,j)=',/,(6e12.4))") bz(:,j)
!     enddo ! j=-1,nlatp2

      end subroutine magdyn
!-----------------------------------------------------------------------
      subroutine dynpotmag
      use fields_module,only: dynpot,phim3d
      use mpi_module,only: mytid
!
! Transform electric potential to geomagnetic grid. 
!
! Electric potential dynpot on geographic grid was defined from source
!   history (see sub rdf4d in nchist.F) before the first step, then by 
!   dynamo module if input flat dynamo=1. Dynpot is dimensioned in 
!   fields_module (fields.F).
!
! Local:
      integer :: k,j,i
      real :: fmin,fmax
!
!     do j=1,nlat
!       do k=1,nlevp1
!         write(6,"('dynpotmag: j=',i2,' k=',i2,' dynpot(:,j,k)=',/,
!    |      (6e12.4))") j,k,dynpot(:,j,k)
!       enddo
!     enddo
!
! Sub geo2mag input dynpot is transformed to magnetic grid in phim3d output.
!     subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlonp1_mag,
!    |  nlon_mag,nlat_mag,lat)
!
!     do k=1,nlevp1
!       do j=1,nmlat
!         call geo2mag(phim3d(1,j,k),dynpot(1,0,k),ig,jg,wt,nlonp1,
!    |      nmlonp1,nmlon,nmlat,j)
!       enddo ! j=1,nmlat
!     enddo ! k=1,nlevp1
!
! Fields.F:
!   dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
!   phim3d(nmlonp1,nmlat,-2:nlevp1), ! 3d electric potential magnetic
! Magfield.F:
!   ig(nmlonp1,nmlat),   ! geog lon index for each geomag point
!   jg(nmlonp1,nmlat),   ! geog lat index for each geomag point
! The "-1" in the j dimension of dynpot is because the dimension is
!   0:nlatp1. Previously, when grdint or geo2mag was called for this
!   transformation, dynpot(1,0,k) was passed as actual arg, and the 
!   dummy arg in the sub was fgeo(nlonp1_geo,*).
!
      do k=1,nlevp1
        do j=1,nmlat
          do i=1,nmlon
            phim3d(i,j,k) =
     |        dynpot(ig(i,j)  ,jg(i,j)  -1,k)*wt(1,i,j)+
     |        dynpot(ig(i,j)+1,jg(i,j)  -1,k)*wt(2,i,j)+
     |        dynpot(ig(i,j)+1,jg(i,j)+1-1,k)*wt(3,i,j)+
     |        dynpot(ig(i,j)  ,jg(i,j)+1-1,k)*wt(4,i,j)
          enddo ! i=1,nmlon
          phim3d(nmlonp1,j,k) = phim3d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=1,nlevp1
!
! Periodic point:
      do k=1,nlevp1
        do j=1,nmlat
          phim3d(nmlonp1,j,k) = phim3d(1,j,k)
        enddo ! j=1,nmlat
      enddo ! k=1,nlevp1
      end subroutine dynpotmag
!-----------------------------------------------------------------------
      subroutine sunloc(iday,secs)
      use cons_module,only: pi,dlamda,dphi
!
! Calculate sun's longitude in dipole coordinates, defining sunlons(nlat)
! (sunlons is module data above). This is called every timestep from advance.
!
! Args:
      integer,intent(in) :: iday ! day of year
      real,intent(in) :: secs    ! ut in seconds
!
! Local:
      integer :: j,isun,jsun
      real :: glats,glons,pisun,pjsun,sndlons,csdlons
!
! Sun's geographic coordinates:
      glats=asin(.398749*sin(2.*PI*(iday-80)/365.))
      glons=pi*(1.-2.*secs/86400.)
!
      pisun = (glons+pi)/dlamda+1.
      pjsun = (glats+.5*(pi-dphi))/dphi+1.
      isun = int(pisun)
      jsun = int(pjsun)
      pisun = pisun-float(isun)
      pjsun = pjsun-float(jsun)
      sndlons = (1.-pisun)*(1.-pjsun)*sin(rlonm(isun+2,jsun))+
     |           pisun*(1.-pjsun)*sin(rlonm(isun+3,jsun))+
     |           pisun*pjsun*sin(rlonm(isun+3,jsun+1))+
     |           (1.-pisun)*pjsun*sin(rlonm(isun+2,jsun+1))
      csdlons = (1.-pisun)*(1.-pjsun)*cos(rlonm(isun+2,jsun))+
     |           pisun*(1.-pjsun)*cos(rlonm(isun+3,jsun))+
     |           pisun*pjsun*cos(rlonm(isun+3,jsun+1))+
     |           (1.-pisun)*pjsun*cos(rlonm(isun+2,jsun+1))
      sunlons(1) = atan2(sndlons,csdlons)
      do j = 2,nlat
        sunlons(j) = sunlons(1)
      enddo
!     write(6,"('sunloc: sunlons=',/,(6e12.4))") sunlons
      end subroutine sunloc
!-----------------------------------------------------------------------
!     subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlonp1_mag,
!    |  nlon_mag,nlat_mag,lat)
!
! Transform field fgeo on geographic grid to geomagnetic grid using
!   indices long,latg and weights wght. Return field fmag on magnetic grid.
! This routine is similiar to sub geo2mag in dynamo.F.
!
! Args:
!     integer,intent(in) :: nlonp1_geo,nlonp1_mag,nlon_mag,nlat_mag,lat
!     integer,dimension(nlonp1_mag,nlat_mag),intent(in) :: long,latg
!     real,intent(in) :: fgeo(nlonp1_geo,*),wght(4,nlonp1_mag,nlat_mag)
!     real,intent(out) :: fmag(nlonp1_mag,*)
!     integer,intent(in) :: iprint
!
! Local:
!     integer :: i
!
!     do i=1,nlon_mag
!       fmag(i,1) =
!    |    fgeo(long(i,lat)  ,latg(i,lat)  )*wght(1,i,lat)+
!    |    fgeo(long(i,lat)+1,latg(i,lat)  )*wght(2,i,lat)+
!    |    fgeo(long(i,lat)+1,latg(i,lat)+1)*wght(3,i,lat)+
!    |    fgeo(long(i,lat)  ,latg(i,lat)+1)*wght(4,i,lat)
!       if (iprint > 0) write(6,"('geo2mag: i=',i3,' lat=',i3,' long=',
!    |    i3,' latg=',i3,' wght=',4e12.4,' fgeo=',e12.4,' fmag=',
!    |    e12.4)") i,lat,long(i,lat),latg(i,lat),wght(:,i,lat),
!    |    fgeo(long(i,lat),latg(i,lat)),fmag(i,1)
!     enddo
!     end subroutine geo2mag
!-----------------------------------------------------------------------
      end module magfield_module
