!
      subroutine readsource(ier)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Read source history.
!
      use params_module,only: nlevp1
      use input_module,only: source,source_start,output,pristart,
     |  f107,f107a,power,noutput
      use hist_module,only: nsource,nhist,ioutfile,ncid
      use fields_module,only: tlbc,ulbc,vlbc,tlbc_nm,ulbc_nm,vlbc_nm,
     |  he,he_nm,field_initval,ar,ar_nm,n2d,n2d_nm,gzigm1,gzigm2,gnsrhs,
     |  zigm1,zigm2,itp
      use ar_module,only: ar_glbm
      use mpi_module,only: lat0,lat1,lon0,lon1,mlat0,mlat1,mlon0,mlon1,
     |  mp_periodic_f4d,mp_periodic_f2d,mp_bndlons_f2d,mp_bndlats_f2d
      use pdynamo_module,only: define_phim3d,geo2mag_2d
      use lbc,only: pshelb
      use my_esmf,only: geo_adota1,geo_adota2,mag_adota1,mag_adota2
      implicit none
!
! Arg:
      integer,intent(out) :: ier
!
! Local:
      integer :: nth,i,k
      real :: flbc(lon0:lon1,lat0:lat1,9)
      real :: flbc_2(lon0-2:lon1+2,lat0-2:lat1+2,9)
!
      ier = 0
!
! If source file was provided, open read-only, and close afterwards.
! If source file not provided, open first output file, and leave it
!   open for possible appending later by output_hist.
!
! Source was provided -- read source history from source file:
      if (nsource==1) then 
        call rdsource(source,source_start,ncid,nth,.false.)
        if (ncid==0) call shutdown('Source history not found')
        nhist  = 0         ! no histories on output file
        ioutfile = 0       ! no output file currently in use
!
! Source file was not provided -- search output files:
      else
        do i = 1,noutput
          call rdsource(output(i),pristart(:,1),ncid,nth,.false.)
          if (ncid > 0) then   ! history was found (not found if ncid==0)
            nhist  = nth       ! number of histories on current output file
            ioutfile = i       ! current output file name is output(ioutfile)
            exit               ! exit the loop if history is found
          endif
        enddo
        if (ncid == 0) then    ! history not found on any output file
          write(6,"('>>> rdsource: starting history at mtime ',4i4,
     |      ' was not found on OUTPUT files.')") pristart(:,1)
          call shutdown('Starting history not found on output files')
        endif
      endif
!
!     if (ncid==0) then
!       ier = 1
!       return
!     endif
!
! Do mpi periodic points exchange for f4d(:) 
! Moved here from sub nc_rdhist because mpi necessary when f4d data
!  is allocated only for task-local subdomain block.
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
      call mp_periodic_f4d(itp)
!
! Periodic points for t,u,v lbc
! These are dimensioned (lond0:lond1,latd0:latd1) (see fields.F):
! real :: flbc(lon0:lon1,lat0:lat1,6)
!
      flbc_2(:,:,1) = tlbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,2) = ulbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,3) = vlbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,4) = tlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,5) = ulbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,6) = vlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,7) = gzigm1(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,8) = gzigm2(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,9) = gnsrhs(lon0-2:lon1+2,lat0-2:lat1+2)

      call mp_bndlons_f2d(flbc_2,lon0,lon1,lat0,lat1,9)

      tlbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,1)
      ulbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,2)
      vlbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,3)
      tlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,4)
      ulbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,5)
      vlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,6)
      gzigm1(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,7)
      gzigm2(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,8)
      gnsrhs(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,9)
!    
      flbc_2(:,:,1) = tlbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,2) = ulbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,3) = vlbc(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,4) = tlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,5) = ulbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,6) = vlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,7) = gzigm1(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,8) = gzigm2(lon0-2:lon1+2,lat0-2:lat1+2)
      flbc_2(:,:,9) = gnsrhs(lon0-2:lon1+2,lat0-2:lat1+2)

      call mp_bndlats_f2d(flbc_2,lon0,lon1,lat0,lat1,9)

      tlbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,1)
      ulbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,2)
      vlbc(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,3)
      tlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,4)
      ulbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,5)
      vlbc_nm(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,6)
      gzigm1(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,7)
      gzigm2(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,8)
      gnsrhs(lon0-2:lon1+2,lat0-2:lat1+2) = flbc_2(:,:,9)
!
      flbc(:,:,1) = tlbc(lon0:lon1,lat0:lat1)
      flbc(:,:,2) = ulbc(lon0:lon1,lat0:lat1)
      flbc(:,:,3) = vlbc(lon0:lon1,lat0:lat1)
      flbc(:,:,4) = tlbc_nm(lon0:lon1,lat0:lat1)
      flbc(:,:,5) = ulbc_nm(lon0:lon1,lat0:lat1)
      flbc(:,:,6) = vlbc_nm(lon0:lon1,lat0:lat1)
      flbc(:,:,7) = gzigm1(lon0:lon1,lat0:lat1)
      flbc(:,:,8) = gzigm2(lon0:lon1,lat0:lat1)
      flbc(:,:,9) = gnsrhs(lon0:lon1,lat0:lat1)

      call mp_periodic_f2d(flbc,lon0,lon1,lat0,lat1,9)

      tlbc(lon0:lon1,lat0:lat1) = flbc(:,:,1)
      ulbc(lon0:lon1,lat0:lat1) = flbc(:,:,2)
      vlbc(lon0:lon1,lat0:lat1) = flbc(:,:,3)
      tlbc_nm(lon0:lon1,lat0:lat1) = flbc(:,:,4)
      ulbc_nm(lon0:lon1,lat0:lat1) = flbc(:,:,5)
      vlbc_nm(lon0:lon1,lat0:lat1) = flbc(:,:,6)
      gzigm1(lon0:lon1,lat0:lat1) = flbc(:,:,7)
      gzigm2(lon0:lon1,lat0:lat1) = flbc(:,:,8)
      gnsrhs(lon0:lon1,lat0:lat1) = flbc(:,:,9)
!
      call define_phim3d(itp)
!
! Get field-line integrated Pedersen/Hall conductances on mag grids,
! used to send conductances (mag) to magnetosphere before the main loop
      call geo2mag_2d(gzigm1(lon0:lon1,lat0:lat1),
     |  zigm1(mlon0:mlon1,mlat0:mlat1),geo_adota1,mag_adota1,'ZIGM1')
      call geo2mag_2d(gzigm2(lon0:lon1,lat0:lat1),
     |  zigm2(mlon0:mlon1,mlat0:mlat1),geo_adota2,mag_adota2,'ZIGM2')
!
! If n2d_nm was not read from the source history, then set it to n2d
!
      if (all(n2d_nm==field_initval)) then
        write(6,"('readsource: N2D_NM apparently not read from source',
     |    ' history. Will init N2D_NM to N2D')")
        n2d_nm = n2d
      endif
!
! If argon was not read from the source history, then set it to
! column global means (see ar_glbm(nlevp1) in comp_ar.F).
!
      if (all(ar==field_initval)) then
        write(6,"('readsource: Argon apparently not read from source',
     |    ' history. Will init AR and AR_NM to ar_glb')")
        do k=1,nlevp1
          ar   (k,:,:,:) = ar_glbm(k)
          ar_nm(k,:,:,:) = ar_glbm(k)
        enddo
      endif
!
! If helium was not read from history, then set it to pshelb globally.
! (pshelb is a constant set in lbc.F). Helium is a 4d prognostic, which
! was initialized to field_initval (probably zero) by sub init_4d (fields.F)
!
      if (all(he==field_initval)) then
        write(6,"('readsource: Helium apparently not read from source',
     |    ' history. Will init HE and HE_NM to pshelb=',e12.4)") pshelb
        he = pshelb    ! whole-array op
        he_nm = pshelb ! whole-array op
      endif 
!
      end subroutine readsource
!-------------------------------------------------------------------
      subroutine rdsource(filepath,mtime,ncid,nthist,reopen_append)
!
! Acquire and read source history:
!
! On input:
!   filepath = mss path to history file
!   mtime(4) = model time of requested source history
!   reopen_append: if true, reopen the file for later writing after
!                  reading the history.
!
! On output:
!   ncid   = file id of history file
!   nthist = source history is nth history on the file 
!
      use nchist_module,only: nc_open,nc_close,nc_rdhist
      implicit none
! 
! Args:
      character(len=*),intent(in) :: filepath
      integer,intent(in) :: mtime(4)
      integer,intent(out) :: nthist,ncid
      logical,intent(in) :: reopen_append
!
! Local:
      integer :: ier              ! error flag
      character(len=120) :: diskfile
!
! Acquire source file:
      diskfile = ' '
      call getfile(filepath,diskfile)
      write(6,"('Acquired source history file ',a,
     |  /,'  (disk file is ',a,')')") trim(filepath),trim(diskfile)
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
        call nc_close(ncid)
        ncid = 0
        return
      endif
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
