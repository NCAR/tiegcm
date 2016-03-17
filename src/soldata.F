!
      module soldata_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use input_module,only: see_ncfile
      implicit none
#ifdef SUN
#include <netcdf3.inc>
#else
#include <netcdf.inc>
#endif
      integer :: 
     |  nstruct,    ! # of structure_elements (unlimited dim)
     |  ndate,      ! # dates
     |  nwave       ! # bins
      integer,allocatable :: date(:)  ! date(ndate)
      real,allocatable :: 
     |  uttime(:,:),   ! UT time in hour since UT=0 of the day (uttime(ndate,nstruct))
     |  yfrac(:,:),    ! UT time in fractional year (yfrac(ndate,nstruct))
     |  wave1(:,:),    ! lower boundary of each bin (wave1(nwave,nstruct))
     |  wave2(:,:),    ! upper boundary of each bin (wave2(nwave,nstruct))
     |  spflux(:,:,:)  ! flux data (nwave,ndate,nstruct)
!
! Current SEE flux, interpolated to current model time:
      real,allocatable :: soldata(:) ! soldata(nwave) 
      real,allocatable :: sflux1(:),sflux2(:)
      contains
!-----------------------------------------------------------------------
      subroutine rd_soldata
!
! Read SEE flux nc file. This is called once per run from init.
!
      integer :: i,ncid,istat,ndims,nvars,ngatts,id_unlim,len,itype,
     |  nuttime,nyfrac,nwave1,nwave2,ndateflux,nwaveflux
      integer ::
     |  idim_date, idv_date,        ! id's for date dimension and var
     |  idim_uttime, idv_uttime,    ! id's for uttime dimension and var
     |  idim_yfrac, idv_yfrac,      ! id's for yfrac dimension and var
     |  idim_wave1, idv_wave1,      ! id's for wave1 dimension and var
     |  idim_wave2, idv_wave2,      ! id's for wave2 dimension and var
     |  idim1_spflux, idim2_spflux, ! dim id's for spflux
     |  idv_spflux
      character(len=120) :: diskfile
      character(len=120) :: attname,atttext
!
      write(6,"(/,72('-'),/,'Read solar data:')")
      diskfile = ' '
      call getfile(see_ncfile,diskfile)
      write(6,"('Acquired SEE solar data file ',a,
     |  /,'  (disk file is ',a,')')") trim(see_ncfile),trim(diskfile) 
!
! Open netcdf dataset:
      istat = nf_open(diskfile,NF_NOWRITE,ncid) 
      if (istat /= NF_NOERR) then
        write(6,"('>>> rd_soldata: error opening see_ncfile ',a)")
     |    trim(diskfile)
        call handle_ncerr(istat,'error opening see_ncfile')
      else
        write(6,"('Opened netCDF file ',a,' istat=',i3,' ncid=',
     |    i3)") trim(diskfile),istat,ncid
      endif
!
! Get number of dims, vars, atts. Print global attributes to stdout.
      istat = nf_inq(ncid,ndims,nvars,ngatts,id_unlim)
      write(6,"(/,'Global attributes from file ',a,':')") 
     |  trim(diskfile)
      do i=1,ngatts
        attname = ' ' ; atttext = ' '
        istat = nf_inq_attname(ncid,NF_GLOBAL,i,attname)
        istat = nf_inq_att(ncid,NF_GLOBAL,attname,itype,len)
        if (itype==NF_CHAR) then
          istat = nf_get_att_text(ncid,NF_GLOBAL,attname,atttext)
          write(6,"('  ',a,': ',a)") trim(attname),trim(atttext)
        endif
      enddo
!
! Get length of unlimited dimension (nstruct):
      istat = nf_inq_dimlen(ncid,id_unlim,nstruct)
!
! Get ndate dimension (dim1_DATE):
      istat = nf_inq_dimid(ncid,'dim1_DATE',idim_date)
      istat = nf_inq_dimlen(ncid,idim_date,ndate)
!
! Get nuttime dimension (dim1_UTTIME):
      istat = nf_inq_dimid(ncid,'dim1_UTTIME',idim_uttime)
      istat = nf_inq_dimlen(ncid,idim_uttime,nuttime)
!
! Get nyfrac dimension (dim1_YFRAC):
      istat = nf_inq_dimid(ncid,'dim1_YFRAC',idim_yfrac)
      istat = nf_inq_dimlen(ncid,idim_yfrac,nyfrac)
!
! Get nwave1 dimension (dim1_WAVE1):
      istat = nf_inq_dimid(ncid,'dim1_WAVE1',idim_wave1)
      istat = nf_inq_dimlen(ncid,idim_wave1,nwave1)
!
! Get nwave2 dimension (dim1_WAVE2):
      istat = nf_inq_dimid(ncid,'dim1_WAVE2',idim_wave2)
      istat = nf_inq_dimlen(ncid,idim_wave2,nwave2)
!
! Get dimensions for spflux:
      istat = nf_inq_dimid(ncid,'dim1_SP_FLUX',idim1_spflux)
      istat = nf_inq_dimlen(ncid,idim1_spflux,nwaveflux)
      istat = nf_inq_dimid(ncid,'dim2_SP_FLUX',idim2_spflux)
      istat = nf_inq_dimlen(ncid,idim2_spflux,ndateflux)
!
!     write(6,"(/,'rd_soldata: nstruct=',i3,' ndate=',i3,' nwave1=',i3,
!    |  ' nwave2=',i3,' nwaveflux=',i3,' ndateflux=',i3)")
!    |  nstruct,ndate,nwave1,nwave2,nwaveflux,ndateflux
!
! Check dimensions:
      if (ndate /= ndateflux) then
        write(6,"('>>> rd_soldata: ndate /= ndateflux: ndate=',i4,
     |    ' ndateflux=',i4)") ndate,ndateflux
        call shutdown('rd_soldata')
      endif
      if (nwave1 /= nwave2) then
        write(6,"('>>> rd_soldata: nwave1 /= nwave2: nwave1=',i4,
     |    ' nwave2=',i4)") nwave1,nwave2
        call shutdown('rd_soldata')
      endif
      nwave = nwave1
      if (nwaveflux /= nwave) then
        write(6,"('>>> rd_soldata: nwaveflux /= nwave: nwaveflux=',i4,
     |    ' nwave=',i4)") nwaveflux,nwave
        call shutdown('rd_soldata')
      endif
      if (nuttime /= ndate) then
        write(6,"('>>> rd_soldata: nuttime /= ndate: nuttime=',i4,
     |    ' ndate=',i4)") nuttime,ndate
        call shutdown('rd_soldata')
      endif
      if (nyfrac /= ndate) then
        write(6,"('>>> rd_soldata: nyfrac /= ndate: nyfrac=',i4,
     |    ' ndate=',i4)") nyfrac,ndate
        call shutdown('rd_soldata')
      endif
!
! Allocate and read date variable:
      if (allocated(date)) deallocate(date)
      allocate(date(ndate),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' date var: ndate=',i4)") ndate
      istat = nf_inq_varid(ncid,'DATE',idv_date)
      istat = nf_get_var_int(ncid,idv_date,date)
      write(6,"('DATE min,max=',2i8)") minval(date),maxval(date)
!
! Allocate and read uttime:
      if (allocated(uttime)) deallocate(uttime)
      allocate(uttime(ndate,nstruct),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' uttime var: ndate=',i4,' nstruct=',i4)") ndate,nstruct
      istat = nf_inq_varid(ncid,'UTTIME',idv_uttime)
      istat = nf_get_var_double(ncid,idv_uttime,uttime)
      write(6,"('UTTIME min,max=',2e12.4)") minval(uttime),
     |            maxval(uttime)
!
! Allocate and read yfrac:
      if (allocated(yfrac)) deallocate(yfrac)
      allocate(yfrac(ndate,nstruct),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' yfrac var: ndate=',i4,' nstruct=',i4)") ndate,nstruct
      istat = nf_inq_varid(ncid,'YFRAC',idv_yfrac)
      istat = nf_get_var_double(ncid,idv_yfrac,yfrac)
      write(6,"('YFRAC min,max=',2e12.4)") minval(yfrac),maxval(yfrac)
!
! Allocate and read wave1:
      nwave = nwave1
      if (allocated(wave1)) deallocate(wave1)
      allocate(wave1(nwave,nstruct),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' wave1 var: nwave=',i4,' nstruct=',i4)") nwave,nstruct
      istat = nf_inq_varid(ncid,'WAVE1',idv_wave1)
      istat = nf_get_var_double(ncid,idv_wave1,wave1)
      write(6,"('WAVE1 min,max=',2e12.4)") minval(wave1),maxval(wave1)
!
! Allocate and read wave2:
      if (allocated(wave2)) deallocate(wave2)
      allocate(wave2(nwave,nstruct),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' wave2 var: nwave=',i4,' nstruct=',i4)") nwave,nstruct
      istat = nf_inq_varid(ncid,'WAVE2',idv_wave2)
      istat = nf_get_var_double(ncid,idv_wave2,wave2)
      write(6,"('WAVE2 min,max=',2e12.4)") minval(wave2),maxval(wave2)
!
! Allocate and read spflux:
      if (allocated(spflux)) deallocate(spflux)
      allocate(spflux(nwave,ndateflux,nstruct),stat=istat)
      if (istat /= 0) write(6,"('>>> rd_soldata: error allocating ',
     |  ' spflux var: nwave=',i4,' ndateflux=',i4,' nstruct=',i3)") 
     |  nwave,ndateflux,nstruct
      istat = nf_inq_varid(ncid,'SP_FLUX',idv_spflux)
      istat = nf_get_var_double(ncid,idv_spflux,spflux)
      write(6,"('SP_FLUX min,max=',2e12.4)") 
     |  minval(spflux),maxval(spflux)
!
! Allocate array for current solar flux:
      if (allocated(soldata)) deallocate(soldata)
      allocate(soldata(nwave),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> rd_soldata: error allocating ',
     |    ' soldata var: nwave=',i4)") nwave
      endif
      if (allocated(sflux1)) deallocate(sflux1)
      allocate(sflux1(nwave),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> rd_soldata: error allocating ',
     |    ' sflux1 var: nwave=',i4)") nwave
      endif
      if (allocated(sflux2)) deallocate(sflux2)
      allocate(sflux2(nwave),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> rd_soldata: error allocating ',
     |    ' sflux2 var: nwave=',i4)") nwave
      endif
!
      write(6,"('Finished reading solar flux data',/,72('-'))")
      end subroutine rd_soldata
!-------------------------------------------------------------------
      subroutine get_soldata(iyear,iday,iutsec,istep)
!
! Args:
      integer,intent(in) :: iyear,iday,iutsec,istep
!
! Local:
      integer :: index,index1,index2,i
      real :: iuttime,iyfrac
      integer,parameter :: istruct = 1
!
! External:
      integer,external::real_bsearch
!
!     write(6,"('get_soldata: iyear=',i5,' iday=',i4,' iutsec=',
!    |  i8)") iyear,iday,iutsec
!
! calculate the fractional year value for the model time
      iuttime= iutsec/3600.
      if (mod(iyear,4) == 0) then
         iyfrac=iyear+(iday-1+(iuttime/24.))/366.
      else
         iyfrac=iyear+(iday-1+(iuttime/24.))/365.
      endif
!
! Check that requested date is available:
      if (iyfrac < yfrac(1,istruct) .or. 
     |        iyfrac > yfrac(ndate,istruct)) then
        write(6,"(/,'>>> get_soldata: requested year-day is not ',
     |    'available from see_ncfile data file:')")
        write(6,"(4x,' date(1)=',i8,
     |    ' date(ndate)=',i8,/)") date(1),date(ndate)
        write(6,"(72('-'),/)")
        call shutdown('soldata')
      endif
!
! Get index to requested time:

      index=real_bsearch(yfrac(:,istruct),1,ndate,iyfrac)
      if (index == -1) then
         call shutdown('soldata: no SEE data available for the date')
      endif

!      write(12,"('index=',i4,'date(index)=',i8)") index,date(index)
!
! make sure that the solar fluxes of both sides of the bracket are not 
! missing (For now, we found that some spectrum is a zero array)
      index1=index
      do
        if (maxval(spflux(:,index1,istruct)) ==0) then
	   index1=index1-1
	   if (index1 < 1 ) then
             write(6,"('>>> No non zero solar spectrum available ',
     |         'for interpolation')")
             call shutdown('soldata')
           endif
        else
	   exit
        endif
      end do
      index2=index+1
      do
        if (maxval(spflux(:,index2,istruct)) ==0) then
	   index2=index2+1
	   if (index2 >ndate ) then
             write(6,"('>>>No non zero solar spectrum available ',
     |         'for interpolation')")
             call shutdown('soldata')
           endif
        else
	   exit
        endif
      end do
!
      sflux1(:) = spflux(:,index1,istruct) 
      sflux2(:) = spflux(:,index2,istruct)

!       write(6,"('date(index1)=',i8,' uttime1=',i8,' sflux1=',/,
!    |    (6e12.4))") date(index1),uttime(index1),sflux1
!       write(6,"('date(index2)  =',i8,' uttime2=',i8,' sflux2=',/,
!    |    (6e12.4))") date(index2),uttime(index),sflux2

      call finterpb(sflux1,sflux2,soldata,nwave,yfrac(index1,istruct),
     |              yfrac(index2,istruct),iyfrac)
      end subroutine get_soldata
!-------------------------------------------------------------------
      subroutine handle_ncerr(istat,msg)
!
! Handle a netcdf lib error:
!
! (Cannot use-associate this from nchist because of circular dependency
!  problems)
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
      end subroutine handle_ncerr
!-----------------------------------------------------------------------
      end module soldata_module
