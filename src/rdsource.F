!
      subroutine readsource(ier)
!
! Read source history.
!
      use input_module,only: source,tempdir,source_start,output,start,
     |  f107,f107a,power,ctpoten,dynamo
      use hist_module,only: nsource,nhist,ioutfile,ncid,iprint,h
#ifdef MPI
      use mpi_module,only: mp_periodic_f4d,mp_dynpot,lat0,lat1,lon0,
     |  lon1
#endif
      use fields_module,only: itp,poten,un,f4d
      implicit none
!
! Arg:
      integer,intent(out) :: ier
!
! Local:
      integer :: lu,nth,ncid_source,j,i,k
      real :: fmin,fmax,fminj,fmaxj
!
      ier = 0
!
! If source file was provided, open read-only, and close afterwards.
! If source file not provided, open first output file, and leave it
!   open for possible appending later by output_hist.
!
! Source was provided -- read source history from source file:
      if (nsource==1) then 
        call rdsource(source,tempdir,source_start,ncid,nth,
     |    .false.,iprint)
        nhist  = 0         ! no histories on output file
        ioutfile = 0       ! no output file currently in use
!
! Source file was not provided -- search 1st output file:
      else
        call rdsource(output(1),tempdir,start(:,1),ncid,nth,.true.,
     |    iprint)
        nhist  = nth       ! number of histories on current output file
        ioutfile = 1       ! current output file name is output(ioutfile)
      endif
!
      if (ncid==0) then
        ier = 1
        return
      endif

!     do j=lat0,lat1
!       do i=lon0,lon1
!         write(6,"('after rdsource: itp=',i3,' lat=',i3,
!    |      ' (lat0,1=',2i3,') i=',i3,' (lon0,1=',2i3,') un(:,i)=',
!    |      /,(6e12.4))") itp,j,lat0,lat1,i,lon0,lon1,un(:,i,j,itp)

!         write(6,"('after rdsource: itp=',i3,' lat=',i3,
!    |      ' (lat0,1=',2i3,') i=',i3,' (lon0,1=',2i3,') f4d(2)%data=',
!    |      /,(6e12.4))") itp,j,lat0,lat1,i,lon0,lon1,
!    |      f4d(2)%data(:,i,j,itp)
!       enddo ! i=lon0,lon1
!     enddo ! j=lat0,lat1

!
! Do mpi periodic points exchange for f4d(:) 
! Moved here from sub nc_rdhist because mpi necessary when f4d data
!  is allocated only for task-local subdomain block.
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
#ifdef MPI
      call mp_periodic_f4d(itp)
!
! Define dynpot(nlonp1,0:nlatp1,nlevp1) from electric potential
! subdomains, which were read from source history. Do this only
! if dynamo input flag is set:
!
      if (dynamo > 0) then
        call mp_dynpot(itp)
      else
        poten(:,:,:,itp) = 0. ! zero out electric potential
      endif
#else
      call set_periodic_f4d(itp)      
!
! Set_dynpot for non-mpi runs:
      call set_dynpot(itp)
#endif
!
! Update gpi vars:
!     f107  = h%f107d
!     f107a = h%f107a
!     power = h%hpower
!     ctpoten = h%ctpoten
!
! Debug lat min,max:
!     write(6,"('readsource: call fgmnmx after reading netcdf',
!    |  ' source: ixtimep=',i2)") ixtimep
!     do j=1,nlat
!       call fgmnmx(0,'DUM',1,ndisk,j,ixtimep)
!     enddo
!
! Debug 3d min,max:
!     write(6,"('readsource: call fminmaxspv after reading netcdf',
!    |  ' source: ixtimep=',i2)") ixtimep
!     do i=1,ndisk
!       fmin = 1.e36
!       fmax = -1.e36
!       do j=1,nlat
!         call fminmaxspv(fg(1,ndexa(i+1)+1,j,ixtimep),zimxp*zkmxp,
!    |      fminj,fmaxj,spval)
!         if (fminj < fmin) fmin = fminj    
!         if (fmaxj > fmax) fmax = fmaxj    
!       enddo
!       write(6,"('readsource: Read field ',a,' 3d min,max = ',
!    |    2e12.4)") nflds_lab(i),fmin,fmax
!     enddo
      end subroutine readsource
!-------------------------------------------------------------------
      subroutine rdsource(filepath,tempdir,mtime,ncid,nthist,
     |  reopen_append,iprint)
!
! Acquire and read source history:
!
! On input:
!   filepath = mss path to history file
!   tempdir  = path to a temporary directory
!   mtime(3) = model time of requested source history
!   reopen_append: if true, reopen the file for later writing after
!                  reading the history.
!   iprint: if > 0, report to stdout
!
! On output:
!   ncid   = file id of history file
!   nthist = source history is nth history on the file 
!   global history structure h is defined (see nc_rdhist)
!
      use nchist_module,only: nc_open,nc_close,nc_rdhist
      implicit none
! 
! Args:
      character(len=*),intent(in) :: filepath,tempdir
      integer,intent(in) :: mtime(3),iprint
      integer,intent(out) :: nthist,ncid
      logical,intent(in) :: reopen_append
!
! Local:
      integer ::
     |  mday,mhour,mmin, ! model day,hour,minute from header
     |  j,               ! latitude loop index
     |  ier              ! error flag
      real :: dum,rj
      character(len=120) :: diskfile
!
! Acquire source file:
! For systems without ncar mss (MSS==0), SOURCE given in user input file 
!   should be a local file name, e.g., SOURCE='ROBLE.tgcm15r.peqsm41.nc'
!
#if (MSS == 0)
      diskfile = filepath
      call get_diskfile(diskfile,' ') 
#else
!     call mkdiskflnm(filepath,diskfile)
!     call getms(filepath,diskfile,tempdir,' ')
      diskfile = ' '
      call getfile(filepath,diskfile)
      write(6,"('Acquired source history file ',a,
     |  /,'  (disk file is ',a,')')") trim(filepath),trim(diskfile)
#endif
!
! Open existing netcdf file for read-only:
      call nc_open(ncid,diskfile,'OLD','READ')
      if (ncid==0) then
        write(6,"(/,'>>> rdsource: error opening ',a,' as a ',
     |    'netcdf file.')") trim(diskfile)
!       return
        call shutdown('open netcdf source history file')
      endif
!
! Search for and read the source history:
      call nc_rdhist(ncid,diskfile,mtime,nthist,ier)
      if (ier > 0) then
        write(6,"(/,'>>> ERROR return from nc_rdhist reading')")
        write(6,"('  source file ',a,' from ncid=',i8,' mtime=',
     |    3i4)") trim(diskfile),ncid,mtime
        call shutdown('nc_rdhist')
      endif
      call nc_close(ncid)
!
! Reopen file for writing if necessary:
      if (reopen_append) call nc_open(ncid,diskfile,'OLD','WRITE')
      end subroutine rdsource
!-----------------------------------------------------------------------
      subroutine set_periodic_f4d(itx)
!
! Set periodic points for all f4d fields (serial or non-mpi only):
!
      use params_module,only: nlonp4
      use fields_module,only: f4d,nf4d_hist
      integer,intent(in) :: itx
      integer :: n

! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
      do n=1,nf4d_hist
        f4d(n)%data(:,1,:,itx) = f4d(n)%data(:,nlonp4-3,:,itx)
        f4d(n)%data(:,2,:,itx) = f4d(n)%data(:,nlonp4-2,:,itx)
!
        f4d(n)%data(:,nlonp4-1,:,itx) = f4d(n)%data(:,3,:,itx)
        f4d(n)%data(:,nlonp4  ,:,itx) = f4d(n)%data(:,4,:,itx)
      enddo
      end subroutine set_periodic_f4d
!-----------------------------------------------------------------------
      subroutine set_dynpot(itx)
!
! This is called for non-MPI runs only:
! Define dynpot(nlonp1,0:nlatp1,nlevp1)  (3d electric potential geographic)
!   from 4d field poten(levd0:levd1,lond0:lond1,latd0:latd1,2) that has
!   been read from source history.
!
      use params_module,only: nlat,nlonp4,nlevp1
      use fields_module,only: poten,dynpot
      integer,intent(in) :: itx
!
! Local: 
      integer :: k,i,j
!
      do j=1,nlat
        do i=1,nlonp4 
          do k=1,nlevp1
            dynpot(i,j,k) = poten(k,i,j,itx) 
          enddo ! k=1,nlevp1
        enddo ! i=1,nlonp4 
      enddo ! j=1,nlat

      end subroutine set_dynpot
