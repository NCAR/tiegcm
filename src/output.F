!
      module output
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use nchist_module,only: ! nchist_mod.F
     |  nc_open,     ! subroutine to open a netcdf file
     |  nc_define,   ! subroutine to define vars in a new netcdf file
     |  nc_wrhist,   ! subroutine to write history to an open netcdf file
     |  nc_close,    ! subroutine to close a netcdf file
     |  nc_fileinfo, ! make contents string for current file
     |  def_fsech,
     |  datasource_kp ! defined here

      use dispose_module,only: add_dispose
!
      contains
!-----------------------------------------------------------------------
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
     |  wrprim, newseries_prim,
     |  wrsech, newseries_sech
      integer :: icount
!
! External:
      logical,external :: wrhist
!
! Determine if it is time to write a history:
! (istep and modeltime are input, remaining 6 args are output)
!
      time2write = wrhist(istep,modeltime,
     |  wrprim, newseries_prim, iseries_prim,
     |  wrsech, newseries_sech, iseries_sech)
!
! Write primary history:
      if (wrprim) then
        call output_hist(istep,newseries_prim)
      endif
!
! Write secondary history:
      if (wrsech) then
        call output_sechist(newseries_sech)
      endif
      end subroutine outhist
!-----------------------------------------------------------------------
      subroutine output_hist(istep,newseries)
!
! Write a primary history to current output disk file. 
! If newseries==T, open new file *before* writing the history
!   (and save existing file before opening the new one)
!
      use input_module,only:  ! input_mod.f
     |  output,      ! output file paths provided by input
     |  noutput,     ! number of output files provided by input
     |  mxhist_prim, ! max number of histories on a primary file
     |  hpss_path    ! dispose path for hsi script
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
     |  iseries_prim,! current primary time series (read only)
     |  copyhist,    ! subroutine to copy history structures
     |  h,sh
      use init_module,only: iyear,iday
      implicit none
!
! Args:
      integer,intent(in) :: istep
      logical,intent(in) :: 
     |  newseries  ! if true, open new file before writing history
!
! Local:
      character(len=1024) :: diskfile 
      integer :: iprinthist
      integer,save :: ihist_total=0
      logical :: fullfile
      character(len=1024) :: fileinfo
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

! 1/13/16 btf: ioutfile is set by readsource (rdsource.F)
!     if (ihist_total==0) then       ! first history of a run
!       ioutfile = 0                 ! first history goes on new file
!       if (nsource==0) ioutfile = 1 ! first history on first output file
!     endif
!
      fullfile = .false.
!
! Use >= because source history may already have > mxhist_prim histories.
! In that case, a new file will be opened.
!
      if (nhist >= mxhist_prim) fullfile = .true.

!     write(6,"('Enter output_hist: nhist=',i4,' mxhist_prim=',i4,
!    |  ' ioutfile=',i4,' noutput=',i4,' fullfile=',l1)") 
!    |  nhist,mxhist_prim,ioutfile,noutput,fullfile
!
! Open new file if starting a new time series, or file was filled
! at last write. Otherwise, open existing file to receive history.
! If disposing to hpss, open local disk file constructed from OUTPUT 
! path, otherwise open OUTPUT path as provided by namelist input.
!
      if (newseries.or.fullfile) then  ! open new file
        ioutfile = ioutfile+1
        if (ioutfile > noutput) then
          write(6,"('>>> Out of primary output files: noutput=',i4,
     |      ' Last output file was ',a)") noutput,output(noutput)
          call shutdown('out of primary output files')
        endif
        nhist = 0             ! number of histories on current file
        if (newseries) write(6,"(/,'Starting primary history',
     |    ' time series ',i2,' at model time ',4i4)") iseries_prim,
     |    modeltime
        diskfile = trim(output(ioutfile))
        call nc_open(ncid,diskfile,'REPLACE','WRITE')
        if (ncid==0) then
          write(6,"(/,'>>> output_hist: Error opening new primary ',
     |      'output history file ',a)") trim(diskfile)
          call shutdown('opening new primary output history file')
        else
          write(6,"('Opened new primary history file ',a,
     |      ' (ioutfile=',i3,')')") trim(diskfile),ioutfile
        endif
      else                            ! append to existing file
        diskfile = trim(output(ioutfile))
        call nc_open(ncid,diskfile,'OLD','APPEND')
        if (ncid==0) then
          write(6,"(/,'>>> output_hist: Error opening existing primary',
     |      'output history file',a)") trim(diskfile)
          call shutdown('opening existing primary output history file')
        else
!         write(6,"('Will append to primary history file (ioutfile=',i4,
!    |      ' nhist=',i4,' file=',a)") ioutfile,nhist,trim(diskfile)
        endif
      endif ! newseries
!
! Define current history structure:
! If writing source history, copy sh into h instead of calling define_hist,
!   and before calling nc_define. (sub copyhist is in hist.F)
!
      nhist = nhist+1 ! increment number of histories on current file
      if (istep == 0) then ! echo source history at beginning of initial run
        sh%ihist = nhist
!
! Set source history model time, day and year from user's
! namelist read for first history of the intial run.
! This may be different than the original source history.
! (i.e., SOURCE_START may be different than START and START_DAY,
!  and START_YEAR may be different that year on the source history)
!
        sh%modeltime = modeltime
        sh%year    = iyear
        sh%day     = iday
        call copyhist(sh,h) ! copy source history sh into h
        write(6,"('output_hist: copied source history structure sh to ',
     |    'history structure h (nhist=',i3,')')") nhist
        h%tuv_lbc_intop = 0 ! this code always writes new histories
      else
        call define_hist('primary')
      endif
!
! Define new nc file with history h:
      if (newseries.or.fullfile)
     |  call nc_define(ncid)  ! define new primary history file
!
! Write the history:
      call nc_wrhist(ncid,diskfile,iprinthist)  ! write the history
      ihist_total = ihist_total+1
!
! Report to stdout:
      write(6,"('Wrote primary   history ',i3,',',i2,',',i2,
     |  ' to ',a,' (',i4,' of ',i4,')')") modeltime(1:3),
     |  trim(diskfile),ihist_total,nhist_total
!
! Add contents string to the file (as global attribute), and to the
! hpss dispose script:
      fileinfo = ' '
      call nc_fileinfo(ncid,fileinfo)
      if (len_trim(hpss_path) > 0) 
     |  call add_dispose(trim(diskfile),trim(fileinfo))
      call nc_close(ncid)
      end subroutine output_hist
!-----------------------------------------------------------------------
      subroutine output_sechist(newseries)
!
! Write a secondary history to current output disk file. 
! If newseries==T, open new file *before* writing the history
!   (and save existing file before opening the new one)
!
      use input_module,only:  ! input_mod.f
     |  secout,      ! output file paths provided by input
     |  mxhist_sech, ! max number of histories on a primary file
     |  hpss_path    ! dispose path for hsi script
      use hist_module,only:   ! hist_mod.f
     |  iprint,      ! print flag (read only)
     |  isecout,     ! index to current output file (read/write)
     |  ncidsech,    ! netcdf file pointer (read/write)
     |  nsech,       ! number of histories on current file (read/write)
     |  modeltime,   ! current model time (read only)
     |  nsech_total, ! total sech histories to be written (read only)
     |  nstep,       ! total number of steps this run
     |  iseries_sech ! current secondary time series (read only)
      implicit none
!
! Args:
      logical,intent(in) :: 
     |  newseries  ! if true, open new file before writing history
!
! Local:
      character(len=1024) :: diskfile 
      integer :: iprinthist
      integer,save :: ihist_total=0
      logical :: fullfile
      character(len=1024) :: fileinfo
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
! If disposing to hpss, open local disk file constructed from SECOUT path, 
! otherwise open SECOUT path as provided by namelist input:
!
      if (newseries.or.fullfile) then       ! open new file
        isecout = isecout+1
        nsech = 0             ! number of histories on current file
        if (newseries) write(6,"(/,'Starting secondary history',
     |    ' time series ',i2,' at model time ',4i4)") iseries_sech,
     |    modeltime
        diskfile = trim(secout(isecout))
        call nc_open(ncidsech,diskfile,'REPLACE','WRITE')
        if (ncidsech==0) then
          write(6,"(/,'>>> output_sechist: Error opening new output ',
     |      'file ',a)") trim(diskfile)
          call shutdown('opening new secondary history output file')
        else
          write(6,"('Opened new secondary history file ',a,
     |      ' (isecout=',i4,')')") trim(diskfile),isecout
        endif
      else                                 ! append to existing file
        diskfile = trim(secout(isecout))
        call nc_open(ncidsech,diskfile,'OLD','APPEND')
        if (ncidsech==0) then
          write(6,"(/,'>>> output_sechist: Error opening existing ',
     |      'secondary history output file ',a)") trim(diskfile)
          call shutdown('open existing secondary history output file')
        endif
      endif ! newseries
!
! Define current history structure:
      nsech = nsech+1 ! increment number of histories on current file
      call define_hist('secondary')

!     write(6,"('output_sechist after define_hist: ncidsech=',i3,
!    |  ' nsech=',i2,' newseries=',l1,' fullfile=',l1)") 
!    |  ncidsech,nsech,newseries,fullfile

      if (newseries.or.fullfile) 
     |  call nc_define(ncidsech)  ! define new secondary history file
!
! First history of first secondary file does not define diags because
! addfld has not been called yet. Therefore, if this is the 2nd history
! of the first file, define diags now:
!
      if (isecout==1.and.nsech==2) call def_fsech(ncidsech)
!
! Write the history:
      call nc_wrhist(ncidsech,diskfile,iprinthist)  ! write the history
      ihist_total = ihist_total+1
!
! Report to stdout:
      write(6,"('Wrote secondary history ',i3,',',i2,',',i2,
     |  ' to ',a,' (',i4,' of ',i4,')')") modeltime(1:3),
     |  trim(diskfile),ihist_total,nsech_total
!
! Add contents string to the file (as global attribute), and to the
! hpss dispose script:
      fileinfo = ' '
      call nc_fileinfo(ncidsech,fileinfo)
      if (len_trim(hpss_path) > 0) 
     |  call add_dispose(trim(diskfile),trim(fileinfo))
      call nc_close(ncidsech)
      end subroutine output_sechist
!-----------------------------------------------------------------------
      subroutine define_hist(hist_type)
!
! Define history structure h in hist_module.F, prior to 
!   writing history file.
! On input, hist_type is either "prim" or "sech" for primary or
!   secondary histories.
!
      use params_module,only: tgcm_name,tgcm_version,spval,nlat,
     |  nlon,nlevp1,zibot,zitop,zmbot,zmtop
      use hist_module,only: sh,h,hist_initype,nhist,nsech,modeltime,
     |  isecout,ioutfile,nsource,nsecsource,iseries_prim,iseries_sech
      use input_module,only: output,secout,step,power,ctpoten,kp,
     |  bximf,byimf,bzimf,swvel,swden,al,f107,f107a,tide,tide2,
     |  colfac,source,hist,sechist,calendar_advance,secsource,
     |  gswm_mi_di_ncfile,gswm_mi_sdi_ncfile,gswm_nm_di_ncfile,
     |  gswm_nm_sdi_ncfile,dynamo,source_start,start,
     |  gpi_ncfile,see_ncfile,ncep_ncfile,label,imf_ncfile,
     |  saber_ncfile,tidi_ncfile,joulefac,gpi_ncfile,bgrddata_ncfile,
     |  ctmt_ncfile
      use mpi_module,only: ntask

      use init_module,only: logname,host,system,iyear,iday,istep,
     |  start_mtime
      use fields_module,only: nf4d_hist,f4d,fsechist
      use cons_module,only: p0,grav
      use aurora_module,only: e1,e2,h1,h2,alfac,ec,alfad,ed
      use gpi_module,only: kp_gpi=>fkp
      implicit none
!
! Args:
      character(len=*),intent(in) :: hist_type ! primary or secondary
!
! Local:
      integer :: i,ier,imo,ida,sourceyear,sourceday,sourcemtime(3)
      character(len=8) :: write_date,write_time
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
      h%hist_type  = hist_type
      if (nsource > 0) then ! initial run
        h%run_type = 'initial'
        h%source_file = source
        h%source_mtime = source_start
      else ! continuation run
        h%run_type = 'continuation'
        h%source_file = output(1)
        h%source_mtime = start(:,1)
      endif
!
! Define primary or secondary history structure (hist_mod.f):
      select case(hist_type)
        case ('primary')
          h%ihist = nhist
          h%delhmins = mtime_to_mins(hist(:,iseries_prim))

          h%output_file = output(ioutfile)
          call expand_path(h%output_file)

          h%nflds = nf4d_hist ! is parameter for now
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
        case ('secondary')
          h%ihist = nsech
          h%delhmins = mtime_to_mins(sechist(:,iseries_sech))
          h%output_file = secout(isecout)
          call expand_path(h%output_file)
        case default
          write(6,"('>>> define_hist: unknown hist_type = ',a)") 
     |      hist_type
      end select
!
      h%label   = label
      h%year    = iyear
      h%day     = iday
      h%calendar_advance = calendar_advance
      call datetime(write_date,write_time)
      h%write_date = write_date//' '//write_time
      h%logname = logname
      h%host    = host
      h%system  = system
      h%model_name    = tgcm_name
      h%model_version = tgcm_version
      h%modeltime = modeltime
      h%time = mtime_to_datestr(h%year,h%modeltime,imo,ida,char80)
      h%iter  = mtime_to_nstep(modeltime(1:3),step)
      h%ut    = float(modeltime(2))+float(modeltime(3))/60.+
     |          float(modeltime(4))/3600.
      h%dtide(:)  = tide2(:)
      h%sdtide(:) = tide(:)
!
! If istep==0, then h%f107, etc were defined from source history (sh),
! otherwise (istep > 0), at least one step has been taken and 
! f107, etc were defined either by user or by getgpi.
!
      if (istep == 0) then
        h%f107d     = sh%f107d
        h%f107a     = sh%f107a
        h%hpower    = sh%hpower
        h%ctpoten   = sh%ctpoten
        h%kp        = sh%kp
        h%bximf     = sh%bximf
        h%byimf     = sh%byimf
        h%bzimf     = sh%bzimf
        h%swvel     = sh%swvel
        h%swden     = sh%swden
        h%al        = sh%al
        h%colfac    = sh%colfac
      else
        h%f107d     = f107
        h%f107a     = f107a
        h%hpower    = power
        h%ctpoten   = ctpoten
!
! Kp can be provided by namelist and/or gpi_ncfile. If both were
! specified, the namelist kp takes precedence (see also advance.F).
!
        if (kp /= spval) then ! write namelist kp to history
          h%kp      = kp
          datasource_kp = 'namelist'
        elseif (len_trim(gpi_ncfile) > 0) then ! write kp from gpi file to history
          h%kp      = kp_gpi
          datasource_kp = 'gpi_ncfile'
        else
          h%kp      = spval ! Kp not used
          datasource_kp = '[none]'
        endif

        h%bximf     = bximf
        h%byimf     = byimf
        h%bzimf     = bzimf
        h%swvel     = swvel
        h%swden     = swden
        h%al        = al
        h%colfac    = colfac
      endif
!
! time-gcm only:
!     h%alfa30    = alfa30
!     h%e30       = e30
!     h%alfad2    = alfa_sp ! was alfad2
!     h%ed2       = e_sp    ! was ed2

      h%e1 = e1
      h%e2 = e2
      h%h1 = h1
      h%h2 = h2
      h%alfac = alfac
      h%alfad = alfad
      h%ec = ec
      h%ed = ed

      h%p0_model  = p0       !  5.e-4 ubars (cons.F)
      h%p0        = p0*1.e-3 !  5.e-7 millibars (hPa)
      h%grav      = grav
      h%step      = step

      h%gpi_ncfile= gpi_ncfile
      if (len_trim(gpi_ncfile) == 0) h%gpi_ncfile = '[none]'
      h%see_ncfile= see_ncfile
      if (len_trim(see_ncfile) == 0) h%see_ncfile = '[none]'
      h%ncep_ncfile= ncep_ncfile
      if (len_trim(ncep_ncfile) == 0) h%ncep_ncfile = '[none]'
      h%imf_ncfile= imf_ncfile
      if (len_trim(imf_ncfile) == 0) h%imf_ncfile = '[none]'

      h%gswm_mi_di_ncfile= gswm_mi_di_ncfile
      if (len_trim(gswm_mi_di_ncfile) == 0) 
     |  h%gswm_mi_di_ncfile = '[none]'

      h%gswm_mi_sdi_ncfile= gswm_mi_sdi_ncfile
      if (len_trim(gswm_mi_sdi_ncfile) == 0) 
     |  h%gswm_mi_sdi_ncfile = '[none]'

      h%gswm_nm_di_ncfile= gswm_nm_di_ncfile
      if (len_trim(gswm_nm_di_ncfile) == 0) 
     |  h%gswm_nm_di_ncfile = '[none]'

      h%gswm_nm_sdi_ncfile= gswm_nm_sdi_ncfile
      if (len_trim(gswm_nm_sdi_ncfile) == 0) 
     |  h%gswm_nm_sdi_ncfile = '[none]'

      h%saber_ncfile= saber_ncfile
      if (len_trim(saber_ncfile) == 0) h%saber_ncfile = '[none]'
      h%tidi_ncfile= tidi_ncfile
      if (len_trim(saber_ncfile) == 0) h%tidi_ncfile = '[none]'
      h%bgrddata_ncfile= bgrddata_ncfile
      if (len_trim(bgrddata_ncfile) == 0) h%bgrddata_ncfile = '[none]'
      h%ctmt_ncfile= ctmt_ncfile
      if (len_trim(ctmt_ncfile) == 0) h%ctmt_ncfile = '[none]'

      h%missing_value = spval
      h%nlat = nlat
      h%nlon = nlon
      h%nlev = nlevp1
      h%zmtop = zmtop
      h%zmbot = zmbot
      h%zitop = zitop
      h%zibot = zibot
      h%dynamo = .false.
      if (dynamo > 0) h%dynamo = .true.
      h%ntask_mpi = ntask
      h%joulefac = joulefac
      end subroutine define_hist
      end module output
