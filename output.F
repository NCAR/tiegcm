!
      subroutine outhist(istep,modeltime)
!
! Determine if its time to write a history, and if so call
! output_hist and/or output_sechist to write primary and/or
! secondary histories as required. This routine is called at 
! every time step from advance (or from tstwrhist).
!
      use hist_module,only: iseries_prim,iseries_sech
      use fields_module,only: itc
      use timing_module
      implicit none
!
! Args:
      integer,intent(in) :: istep,modeltime(4)
!
! Local:
      logical :: time2write,
     |  wrprim,  save_prim, newseries_prim,
     |  wrsech,  save_sech, newseries_sech
      integer :: icount
      real :: elapsed
!
! External:
      logical,external :: wrhist
!
! Determine if it is time to write a history:
! (istep and modeltime are input, remaining 6 args are output)
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, save_prim, newseries_prim, iseries_prim,
     |  wrsech, save_sech, newseries_sech, iseries_sech)
!
! Write primary history:
      if (wrprim) then
        call start_timing(icount,'write primary history')
        call output_hist(save_prim,newseries_prim)
        call end_timing(icount,elapsed)
        elapsed_prim = elapsed_prim+elapsed
      endif
!
! Write secondary history:
      if (wrsech) then
        call start_timing(icount,'write secondary history')
        call output_sechist(save_sech,newseries_sech)
        call end_timing(icount,elapsed)
        elapsed_sech = elapsed_sech+elapsed
      endif
      end subroutine outhist
!-----------------------------------------------------------------------
      subroutine output_hist(svfile,newseries)
!
! Write a primary history to current output disk file. 
! If svfile==T, save the file *after* writing the history.
! If newseries==T, open new file *before* writing the history
!   (and save existing file before opening the new one)
!
      use input_module,only:  ! input_mod.f
     |  output,      ! output file paths provided by input
     |  tempdir,     ! temporary scratch directory for saving
     |  dispose,     ! 0/1 for do/donot save file to ncar mss
     |  mxhist_prim  ! max number of histories on a primary file
      use hist_module,only:   ! hist_mod.f
     |  iprint,      ! print flag (read only)
     |  ioutfile,    ! index to current output file (read/write)
     |  ncid,        ! netcdf file pointer (read/write)
     |  nhist,       ! number of histories on current file (read/write)
     |  modeltime,   ! current model time (read only)
     |  nhist_total, ! total primary histories to be written (read only)
     |  nsource,     ! 0/1 number of source files provided (read only)
     |  nsecsource,  ! 0/1 number of secsource files provided (read only)
     |  nstep,       ! total number of steps this run
     |  iseries_prim ! current primary time series (read only)
      use nchist_module,only: ! nchist_mod.F
     |  nc_open,     ! subroutine to open a netcdf file
     |  nc_define,   ! subroutine to define vars in a new netcdf file
     |  nc_wrhist    ! subroutine to write history to an open netcdf file
      use init_module,only: istep ! current time-step
      implicit none
!
! Args:
      logical,intent(in) :: 
     |  newseries, ! if true, open new file before writing history
     |  svfile     ! if true, save current output file after writing history
!
! Local:
      character(len=80) :: diskfile 
      integer :: iprinthist
      integer,save :: ihist_total=0
      logical :: fullfile
!
! Init:
      diskfile = ' '
      if (iprint > 0) write(6,"(/,72('-'),/,'Primary History Output:',
     |  ' model time (day,hour,min) = ',4i4)") modeltime
      if (ihist_total == 0) then
        iprinthist = 1
      else
        iprinthist = iprint
      endif
      if (ihist_total==0) then       ! first history of a run
        ioutfile = 0                 ! first history goes on new file
        if (nsource==0) ioutfile = 1 ! first history on first output file
      endif
      fullfile = .false.
      if (nhist==mxhist_prim) fullfile = .true.
!
! Open new file if starting a new time series, or file was filled
! at last write:
      if (newseries.or.fullfile) then
        ioutfile = ioutfile+1
        nhist = 0             ! number of histories on current file
        if (newseries) write(6,"(/,'Starting primary history',
     |    ' time series ',i2,' at model time ',4i4)") iseries_prim,
     |    modeltime
        call mkdiskflnm(output(ioutfile),diskfile)
        call nc_open(ncid,diskfile,'REPLACE','WRITE')
      endif ! newseries
!
! Define current history structure:
      nhist = nhist+1 ! increment number of histories on current file
      call define_hist('primary')
      if (newseries.or.fullfile) 
     |  call nc_define(ncid)  ! define new file
!
! Write the history:
      call nc_wrhist(ncid,iprinthist)  ! write the history
      ihist_total = ihist_total+1
!
! Report to stdout:
      call mkdiskflnm(output(ioutfile),diskfile)
      write(6,"('Wrote primary   history ',i3,',',i2,',',i2,
     |  ' to ',a,' (',i3,' of ',i3,')')") modeltime(1:3),
     |  trim(diskfile),ihist_total,nhist_total
!
! Save file if necessary:
      if (svfile.or.nhist==mxhist_prim.or.nstep==0) then
        call savefile(ncid,diskfile,output(ioutfile),tempdir,
     |    dispose,.true.,iprint)
      endif
      end subroutine output_hist
!-----------------------------------------------------------------------
      subroutine output_sechist(svfile,newseries)
!
! Write a secondary history to current output disk file. 
! If svfile==T, save the file *after* writing the history.
! If newseries==T, open new file *before* writing the history
!   (and save existing file before opening the new one)
!
      use input_module,only:  ! input_mod.f
     |  secout,      ! output file paths provided by input
     |  tempdir,     ! temporary scratch directory for saving
     |  dispose,     ! 0/1 for do/donot save file to ncar mss
     |  mxhist_sech  ! max number of histories on a primary file
      use hist_module,only:   ! hist_mod.f
     |  iprint,      ! print flag (read only)
     |  isecout,     ! index to current output file (read/write)
     |  ncidsech,    ! netcdf file pointer (read/write)
     |  nsech,       ! number of histories on current file (read/write)
     |  modeltime,   ! current model time (read only)
     |  nsech_total, ! total sech histories to be written (read only)
     |  nstep,       ! total number of steps this run
     |  iseries_sech ! current secondary time series (read only)
      use nchist_module,only: ! nchist_mod.F
     |  nc_open,     ! subroutine to open a netcdf file
     |  nc_define,   ! subroutine to define vars in a new netcdf file
     |  nc_wrhist    ! subroutine to write history to an open netcdf file
      use init_module,only: istep ! current time-step
      implicit none
!
! Args:
      logical,intent(in) :: 
     |  newseries, ! if true, open new file before writing history
     |  svfile     ! if true, save current output file after writing history
!
! Local:
      character(len=80) :: diskfile 
      integer :: iprinthist
      integer,save :: ihist_total=0
      logical :: fullfile
!
! Init:
      diskfile = ' '
      if (iprint > 0) write(6,"(/,72('-'),/,'Secondary History Output:',
     |  ' model time (day,hour,min) = ',3i4)") modeltime
      if (ihist_total == 0) then
        iprinthist = 1
      else
        iprinthist = iprint
      endif
      if (ihist_total==0) isecout = 0 ! first history of a run
      fullfile = .false.
      if (nsech==mxhist_sech) fullfile = .true.
!
! Open new file if starting a new time series, or file was filled
! at last write:
      if (newseries.or.fullfile) then
        isecout = isecout+1
        nsech = 0             ! number of histories on current file
        if (newseries) write(6,"(/,'Starting secondary history',
     |    ' time series ',i2,' at model time ',4i4)") iseries_sech,
     |    modeltime
        call mkdiskflnm(secout(isecout),diskfile)
        call nc_open(ncidsech,diskfile,'REPLACE','WRITE')
      endif ! newseries
!
! Define current history structure:
      nsech = nsech+1 ! increment number of histories on current file
      call define_hist('secondary')

!     write(6,"('output_sechist after define_hist: ncidsech=',i3,
!    |  ' nsech=',i2,' newseries=',l1,' fullfile=',l1)") 
!    |  ncidsech,nsech,newseries,fullfile

      if (newseries.or.fullfile) 
     |  call nc_define(ncidsech)  ! define new file
!
! Write the history:
      call nc_wrhist(ncidsech,iprinthist)  ! write the history
      ihist_total = ihist_total+1
!
! Report to stdout:
      call mkdiskflnm(secout(isecout),diskfile)
      write(6,"('Wrote secondary history ',i3,',',i2,',',i2,
     |  ' to ',a,' (',i3,' of ',i3,')')") modeltime(1:3),
     |  trim(diskfile),ihist_total,nsech_total
!
! Save file if necessary:
      if (svfile.or.nsech==mxhist_sech) then
        call savefile(ncidsech,diskfile,secout(isecout),tempdir,
     |    dispose,.true.,iprint)
      endif
      end subroutine output_sechist
!-----------------------------------------------------------------------
      subroutine savefile(ncid,diskfile,outpath,tempdir,idispose,
     |  reopen,iprint)
      use input_module,only: msreten
      use nchist_module,only: nc_open,nc_close,nc_fileinfo
!
! Save file diskfile to tempdir and/or mss.
!
      implicit none
!
! Args:
      integer,intent(in) :: idispose,iprint
      logical,intent(in) :: reopen
      integer,intent(inout) :: ncid
      character(len=*),intent(in) :: diskfile,outpath,tempdir
!
! Local:
      integer :: istat,ier,msasync
      character(len=80) :: fileinfo
!
! External:
      integer,external :: ilink,iunlink
!
! Make fileinfo string for mss comment field and add as a
! global file attribute:
!
!     write(6,"('savefile calling nc_fileinfo for diskfile=',a)")
!    |  trim(diskfile)
      call nc_fileinfo(ncid,fileinfo)
      write(6,"('History file info: ',a)") trim(fileinfo)
!
! Close file:
      call nc_close(ncid)
!
! Dispose to mss only if dispose flag is set, otherwise 
! save link on tempdir.
!     subroutine putms(mspath,dskfile,tmpdir,lnkfile,wrpass,msscomment,
!                      msasync)
!
      if (idispose==1) then
        msasync = 1 ! suspicious crash (nf_create) on Cray when msasync==1
!       msasync = 0 ! suspicious crash (nf_create) on Cray when msasync==1
        call putms(outpath,diskfile,tempdir,' ','NCARTGCM',fileinfo,
     |    msreten,msasync)
      else
!
! Do not link to tempdir if tempdir is also cwd:
        if (trim(tempdir) /= '.') then
          istat = iunlink(trim(tempdir)//'/'//trim(diskfile),iprint)
!
! Let ilink report to stdout regardless of iprint:
          istat = ilink(diskfile,trim(tempdir)//'/'//trim(diskfile),1)
          if (istat /= 0) 
     |      write(6,"('>>> savefile: error ',i3,' from ilink.')") istat
          if (iprint > 0)
     |      write(6,"('(File was NOT disposed to mss)')")
        else
          write(6,"('NOT linking file ',a,' to tempdir ',
     |      'because tempdir==cwd')") trim(diskfile)
        endif
      endif
!
! Reopen file for appending:
      if (reopen) call nc_open(ncid,diskfile,'OLD','APPEND')
      end subroutine savefile
!-----------------------------------------------------------------------
      subroutine define_hist(type)
!
! Define history structure h in hist_module.F, prior to 
!   writing history file.
! On input, type is either "prim" or "sech" for primary or
!   secondary histories.
!
      use params_module,only: tgcm_name,tgcm_version,spval,nlat,
     |  nlon,nlevp1,zst,zsb
      use hist_module,only: h,hist_initype,nhist,nsech,modeltime,
     |  isecout,ioutfile,nsource,nsecsource,iseries_prim,iseries_sech
      use input_module,only: output,secout,date,step,power,ctpoten,
     |  byimf,f107,f107a,mag,tide,tide2,colfac,source,hist,sechist,
     |  calendar_advance,
! am_09/02: External from magnetosphere:
     |  secsource
      use init_module,only: rundate,runtime,logname,host,system,
     |  iyear,iday,istep,igetgpi,igetgswmdi,igetgswmsdi,start_mtime,
     |  igetgswmnmidi,igetgswmnmisdi
      use fields_module,only: nf4d,f4d,fsech,fsechmag,fsech2d,
     |  fsechmag2d,fsechmagphr2d
      use cons_module,only: p0
      use aurora_module,only: alfa30,e30,alfa_sp,e_sp
      implicit none
!
! Args:
      character(len=*),intent(in) :: type ! primary or secondary
!
! Local:
      integer :: i,ier,imo,ida,sourceyear,sourceday,sourcemtime(3)
      character(len=8) :: writedate,writetime
      character(len=80) :: char80
!
! External:
      integer,external :: mtime_to_mins,mtime_to_nstep
      real,external :: mtime_to_datestr
!
! Init (note hist_initype deallocates h%fnames):
!
      call hist_initype(h,istep)
!
! Define primary or secondary history structure (hist_mod.f):
      select case(type)
        case ('primary')
          h%ihist = nhist
          h%delhmins = mtime_to_mins(hist(:,iseries_prim))
          h%mss_path = output(ioutfile)
          if (nsource > 0) then ! initial run
            h%mss_source = source
          else ! continuation run
            h%mss_source = output(1)
          endif
          if (nsecsource > 0) then ! initial run
            h%mss_secsource = secsource
          else ! continuation run
            h%mss_secsource = ' not defined '
              write(6,"('>>> define_hist: error allocating ',
     |          'h%mss_secsource for second. history')")
          endif
          h%nflds = nf4d ! is parameter for now
          if (associated(h%fnames)) deallocate(h%fnames)
          if (h%nflds > 0) then
            allocate(h%fnames(h%nflds),stat=ier)
            if (ier /= 0) then
              write(6,"('>>> define_hist: error allocating h%fnames',
     |          ' for primary history: h%nflds=',i3)") h%nflds
            else
              do i=1,h%nflds
                h%fnames(i) = f4d(i)%short_name
              enddo
            endif
          endif
          h%nfgeo = h%nflds
          h%nfmag = 0       ! no magnetic fields on primary history
          h%nfgeo2d = 0     ! no geographic 2d fields on primary history
          h%nfmag2d  = 0    ! no magnetic 2d fields on primary history
          h%nfmagphr = 0    ! no magnetospheric 2d fields on primary history
        case ('secondary')
          if (nsource > 0) then ! initial run
            h%mss_source = source
          else ! continuation run
            h%mss_source = output(1)
          endif
          if (nsecsource > 0) then ! initial run
            h%mss_secsource = secsource
          else ! continuation run
            h%mss_secsource = ' not defined '
              write(6,"('>>> define_hist: error allocating ',
     |          'h%mss_secsource for second. history')")
          endif
          h%ihist = nsech
          h%delhmins = mtime_to_mins(sechist(:,iseries_sech))
          h%mss_path = secout(isecout)
          if (nsource > 0) then ! initial run
            h%mss_source = source
          else ! continuation run
            h%mss_source = output(1)
          endif
!
! h%nfgeo:    # geographic fields
! h%nfmag:    # magnetic fields
! h%nfgeo2d:  # geographic 2d fields           
! h%nfmag2d:  # magnetic 2d fields       
! h%nfmagphr: # magnetospheric fields   
!
          h%nfgeo = count(len_trim(fsech%short_name) > 0)
          h%nfmag = count(len_trim(fsechmag%short_name) > 0)
          h%nfgeo2d  = count(len_trim(fsech2d%short_name) > 0)
          h%nfmag2d  = count(len_trim(fsechmag2d%short_name) > 0)   
          h%nfmagphr = count(len_trim(fsechmagphr2d%short_name) > 0)
          h%nflds = h%nfgeo+h%nfmag+h%nfgeo2d+h%nfmag2d+h%nfmagphr
          if (h%nflds > 0 .and. .not.associated(h%fnames)) then
!           write(6,"('define_hist allocating sech h%fnames with ',
!    |        'h%nflds=',i3)") h%nflds
            allocate(h%fnames(h%nflds),stat=ier)
            if (ier /= 0) then
              write(6,"('>>> define_hist: error allocating h%fnames',
     |          ' for secondary history: h%nflds=',i3)") h%nflds
            else
              if (h%nfgeo > 0) then
                do i=1,h%nfgeo
                  h%fnames(i) = fsech(i)%short_name
                enddo
              endif
              if (h%nfmag > 0) then
                do i=1,h%nfmag
                  h%fnames(h%nfgeo+i) = fsechmag(i)%short_name
                enddo
              endif
              if (h%nfgeo2d > 0) then
                do i=1,h%nfgeo2d
                  h%fnames(h%nfgeo+h%nfmag+i) =fsech2d(i)%short_name
                enddo
              endif
              if (h%nfmag2d > 0) then
                do i=1,h%nfmag2d
                  h%fnames(h%nfgeo+h%nfmag+h%nfgeo2d+i) =
     |                fsechmag2d(i)%short_name
                enddo
              endif
              if (h%nfmagphr > 0) then
                do i=1,h%nfmagphr
                  h%fnames(h%nfgeo+h%nfmag+h%nfgeo2d+h%nfmag2d+i) =
     |                fsechmagphr2d(i)%short_name
                enddo
              endif
            endif
          endif
        case default
          write(6,"('>>> define_hist: unknown type = ',a)") type
      end select
!
      h%year    = iyear
!     h%calday  = iday
      h%day     = iday
      h%calendar_advance = calendar_advance
!
! Rundate and runtime were set by init and are saved in init module.
! H%createdate will be defined by nc_wrhist.
! Writedate is defined here.
! 
      h%rundate = rundate//' '//runtime
      call datetime(writedate,writetime)
      h%writedate = writedate//' '//writetime
      h%logname = logname
      h%host    = host
      h%system  = system
      h%model_name    = tgcm_name
      h%model_version = tgcm_version
      h%type  = type
      h%modeltime = modeltime
!     h%time  = mtime_to_mins(modeltime(1:3))
      h%time = mtime_to_datestr(h%year,h%modeltime,imo,ida,char80)
      h%iter  = mtime_to_nstep(modeltime(1:3),step)
      h%ut    = float(modeltime(2))+float(modeltime(3))/60.+
     |          float(modeltime(4))/3600.
      h%mag(:,:)  = mag(:,:)
      h%dtide(:)  = tide2(:)
      h%sdtide(:) = tide(:)
!
! If istep==0, then h%f107, etc were defined from source history,
! otherwise (istep > 0), at least one step has been taken and 
! f107, etc were defined either by user or by getgpi.
!
      if (istep > 0) then
        h%f107d     = f107
        h%f107a     = f107a
        h%hpower    = power
        h%ctpoten   = ctpoten
      endif
      h%byimf     = byimf
      h%colfac    = colfac
! time-gcm only:
      h%alfa30    = alfa30
      h%e30       = e30
      h%alfad2    = alfa_sp ! was alfad2
      h%ed2       = e_sp    ! was ed2
      h%p0        = p0
      h%step      = step
      h%gpi       = igetgpi
      h%gswmdi    = igetgswmdi
      h%gswmsdi   = igetgswmsdi
      h%gswmnmidi = igetgswmnmidi
      h%gswmnmisdi= igetgswmnmisdi
      h%missing_value = spval
!
      h%nlat = nlat
      h%nlon = nlon
      h%nlev = nlevp1
      h%zptop = zst
      h%zpbot = zsb
      end subroutine define_hist
