      module imf_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use input_module,only: imf_ncfile
      use hist_module,only: modeltime
      implicit none
#include <netcdf.inc>
!
! All are allocated ndata except for f107d and f107a, which
! are allocated ndays.
! 
      integer :: ndata,ndays
      real :: missing
      integer,dimension(:),target,allocatable :: days
      real,dimension(:),target,allocatable :: date,
     |  uttime,yfrac,bx,by,bz,swvel,swden,kp,f107d,f107a

      contains
!-----------------------------------------------------------------------
      subroutine getimf(iyear,iday,iutsec,istep,iprint)
      use input_module,only: bximf,byimf,bzimf,swvel_inp=>swvel,
     |  swden_inp=>swden,f107_inp=>f107,f107a_inp=>f107a,
     |  rd_f107,rd_f107a,rd_bximf,rd_byimf,rd_bzimf,rd_swvel,rd_swden
      use params_module,only: spval
      implicit none
!
! Args:
      integer,intent(in) :: iyear,iday,iutsec,istep,iprint
!
! Local:
      real :: model_uttime,model_yfrac,model_dfrac,ut
      integer :: index,index1,index2,i,iyyyy4
      real :: model_bx,model_by,model_bz,model_swvel,model_swden,
     |  model_kp,model_f107d,model_f107a
      integer(kind=8) :: iyearday,iyyyy,iddd
!
! External:
      integer,external :: real_bsearch,long_bsearch
!
! Read file if first step:
      if (istep==1) call rdimf(istep)
!
! Calculate model ut, fractional year-day, and fractional year:
      model_uttime = iutsec/3600.
      model_dfrac = iyear*1000+real(iday)+model_uttime/24.
      if (mod(iyear,4) == 0) then
         model_yfrac=iyear+(iday-1+(model_uttime/24.))/366.
      else
         model_yfrac=iyear+(iday-1+(model_uttime/24.))/365.
      endif
!
! Calculate fractional year yfrac from data date (yyyyddd.dayfrac):
      do i=1,ndata
        iyearday = int(date(i),8)           ! yyyyddd (long int)
        iyyyy = iyearday/1000               ! yyyy
        iddd = iyearday-iyyyy*1000          ! ddd
        ut = (date(i)-real(iyearday))*24.  ! ut (hrs)
        iyyyy4 = iyyyy
        if (mod(iyyyy4,4) == 0) then
           yfrac(i) = real(iyyyy)+(real(iddd-1)+ut/24.)/366.
        else
           yfrac(i) = real(iyyyy)+(real(iddd-1)+ut/24.)/365.
        endif
!       write(6,"('getimf: i=',i6,' ndata=',i6,' date(i)=',f14.4,
!    |    ' year=',i4,' day=',i3,' ut=',f6.2,' yfrac(i)=',f9.4)") 
!    |    i,ndata,date(i),iyyyy,iddd,ut,yfrac(i)
      enddo
!
! Check that requested date is available:
      if (model_yfrac < yfrac(1) .or.
     |    model_yfrac > yfrac(ndata)) then
        write(6,"(/,'>>> getimf: requested date ',f14.4,' is not ',
     |    'available from imf_ncfile data file:')") model_yfrac
        write(6,"(4x,' date(1)=',f14.4,
     |    ' date(ndata)=',f14.4,/)") date(1),date(ndata)
        call shutdown('getimf')
      endif
!
! Get index to requested time:
      index=real_bsearch(yfrac,1,ndata,model_yfrac)
      if (index == -1) then
         write(6,"('>>> getimf: error from real_bsearch: ',
     |     'could not find model_yfrac=',i8)") model_yfrac
         call shutdown('getimf: IMF data not available')
      endif
      index1=index
      index2=index+1
      if (index2 <= ndata) then
!
! Check for missing data:
! As of 6/2/08, missing data in imf data file is fatal, unless
!   the user has provided namelist read values (constant or
!   time-dependent): 
!
        if ((ismissing(bx(index1)).or.ismissing(bx(index2))).and.
     |       ismissing(rd_bximf)) then
          call printmissing(index1,index2,'BX')
          call shutdown('Missing BX')
        endif
        if ((ismissing(by(index1)).or.ismissing(by(index2))).and.
     |       ismissing(rd_byimf)) then
          call printmissing(index1,index2,'BY')
          call shutdown('Missing BY')
        endif
        if ((ismissing(bz(index1)).or.ismissing(bz(index2))).and.
     |       ismissing(rd_bzimf)) then
          call printmissing(index1,index2,'BZ')
          call shutdown('Missing BZ')
        endif
        if ((ismissing(swvel(index1)).or.ismissing(swvel(index2))).and.
     |       ismissing(rd_swvel)) then
          call printmissing(index1,index2,'SWVEL')
          call shutdown('Missing SWVEL')
        endif
        if ((ismissing(swden(index1)).or.ismissing(swden(index2))).and.
     |       ismissing(rd_swden)) then
          call printmissing(index1,index2,'SWDEN')
          call shutdown('Missing SWDEN')
        endif
!
! Interpolate data to model time:
!
! bx
        call finterpb(bx(index1),bx(index2),model_bx,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
! by
        call finterpb(by(index1),by(index2),model_by,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
! bz
        call finterpb(bz(index1),bz(index2),model_bz,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
! swvel
        call finterpb(swvel(index1),swvel(index2),model_swvel,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
! swden
        call finterpb(swden(index1),swden(index2),model_swden,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
! kp
        call finterpb(kp(index1),kp(index2),model_kp,1,
     |    yfrac(index1),yfrac(index2),model_yfrac)
      endif
!
! Get index to requested day for solar flux:
!     index=real_bsearch(real(days),1,ndays,model_dfrac)
!     if (index == -1) then
!        write(6,"('>>> getimf: error from real_bsearch: ',
!    |     'could not find model_dfrac=',f10.3,' ndays=',i5,
!    |     ' days(1)=',i8,' days(ndays)=',i8)") 
!    |     model_dfrac,ndays,days(1),days(ndays)
!        call shutdown('getimf: IMF data not available')
!     endif
!     index1=index
!     index2=index+1
!     if (index2 <= ndays) then
!
! Interpolate solar flux:
!
! f107d
!       call finterpb(f107d(index1),f107d(index2),model_f107d,1,
!    |    real(days(index1)),real(days(index2)),model_dfrac)
! f107a
!       call finterpb(f107a(index1),f107a(index2),model_f107a,1,
!    |    real(days(index1)),real(days(index2)),model_dfrac)
!     endif
!
! Transfer new interpolated data to model variables:
! (input has already determined that at least one of these
!  parameters was not provided by the user via namelist read)
!
      if (rd_bximf == spval) bximf = model_bx
      if (rd_byimf == spval) byimf = model_by
      if (rd_bzimf == spval) bzimf = model_bz
      if (rd_swvel == spval) swvel_inp = model_swvel
      if (rd_swden == spval) swden_inp = model_swden
!     if (rd_f107  == spval) f107_inp  = model_f107d
!     if (rd_f107a == spval) f107a_inp = model_f107a

!     write(6,"('getimf: data interpolated to modeltime ',3i4,':')")
!    |  modeltime(1:3)
!     write(6,"('  bximf,byimf,bzimf = ',3f8.2,' swvel=',f8.2,
!    |  ' swden=',f8.2,' f107d,f107a=',2f8.2)") bximf,byimf,bzimf,
!    |  swvel_inp,swden_inp,f107_inp,f107a_inp

      end subroutine getimf
!-----------------------------------------------------------------------
      subroutine rdimf(istep)
      use nchist_module,only: handle_ncerr
      use input_module,only: mxlen_filename
!
! Local:
      integer :: istat,ncid,istep
      character(len=240) :: char240
      character(len=mxlen_filename) :: dskfile
      integer :: id_ndata,id_ndays,idv_date,idv_uttime,idv_yfrac,
     |  idv_bx,idv_by,idv_bz,idv_swvel,idv_swden,idv_kp,
     |  idv_f107d,idv_f107a,idv_missing,idv_days
      real :: fmin,fmax
!
      dskfile = ' '
      call getfile(imf_ncfile,dskfile)
!
! Open imf data file:
      istat = nf_open(dskfile,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) then
        write(char240,"('Error opening imf_ncfile ',a)")trim(imf_ncfile)
        call handle_ncerr(istat,char240)
        call shutdown('imf_ncfile')  
      endif
      write(6,"(/,72('-'))")
      write(6,"('RDIMF: read IMF data file:')")
      write(6,"('Opened netcdf imf data file ',a)") trim(dskfile)
!
! Get ndata and ndays dimensions:
      istat = nf_inq_dimid(ncid,'ndata',id_ndata)
      istat = nf_inq_dimlen(ncid,id_ndata,ndata)
!     istat = nf_inq_dimid(ncid,'ndays',id_ndays)
!     istat = nf_inq_dimlen(ncid,id_ndays,ndays)
!     write(6,"('rdimf: ndata=',i5,' ndays=',i5)") ndata,ndays
!
! Allocate data arrays and initialize structure pointers:
      if (istep==1) call alloc_imf
!
! Read vars:
!
! Missing value:
!     istat = nf_inq_varid(ncid,'missing',idv_missing)
!     istat = nf_get_var_double(ncid,idv_missing,missing)

! yyyyddd.frac date:
      istat = nf_inq_varid(ncid,'date',idv_date)
      istat = nf_get_var_double(ncid,idv_date,date)
      write(6,"('rdimf: date min,max=',2f14.4)") 
     |  minval(date),maxval(date)

! integer days:
!     istat = nf_inq_varid(ncid,'days',idv_days)
!     istat = nf_get_var_int(ncid,idv_days,days)
!     write(6,"('rdimf: ndays=',i5,' days min,max=',2i10)") 
!    |  ndays,minval(days),maxval(days)

! uttime:
!     istat = nf_inq_varid(ncid,'uttime',idv_uttime)
!     istat = nf_get_var_double(ncid,idv_uttime,uttime)
!     write(6,"('rdimf: uttime min,max=',2f6.2)") 
!    |  minval(uttime),maxval(uttime)

! yfrac:
!     istat = nf_inq_varid(ncid,'yfrac',idv_yfrac)
!     istat = nf_get_var_double(ncid,idv_yfrac,yfrac)
!     write(6,"('rdimf: yfrac min,max=',2f10.4)") 
!    |  minval(yfrac),maxval(yfrac)

! bx:
      istat = nf_inq_varid(ncid,'bx',idv_bx)
      istat = nf_get_var_double(ncid,idv_bx,bx)
      call fminmaxspv(bx,ndata,fmin,fmax,missing)
      write(6,"('rdimf: bx min,max=',2f8.2)") fmin,fmax

! by:
      istat = nf_inq_varid(ncid,'by',idv_by)
      istat = nf_get_var_double(ncid,idv_by,by)
      call fminmaxspv(by,ndata,fmin,fmax,missing)
      write(6,"('rdimf: by min,max=',2f8.2)") fmin,fmax

! bz:
      istat = nf_inq_varid(ncid,'bz',idv_bz)
      istat = nf_get_var_double(ncid,idv_bz,bz)
      call fminmaxspv(bz,ndata,fmin,fmax,missing)
      write(6,"('rdimf: bz min,max=',2f8.2)") fmin,fmax

! swvel:
      istat = nf_inq_varid(ncid,'swvel',idv_swvel)
      istat = nf_get_var_double(ncid,idv_swvel,swvel)
      call fminmaxspv(swvel,ndata,fmin,fmax,missing)
      write(6,"('rdimf: swvel min,max=',2f8.2)") fmin,fmax

! swden:
      istat = nf_inq_varid(ncid,'swden',idv_swden)
      istat = nf_get_var_double(ncid,idv_swden,swden)
      call fminmaxspv(swden,ndata,fmin,fmax,missing)
      write(6,"('rdimf: swden min,max=',2f8.2)") fmin,fmax

! kp:
!     istat = nf_inq_varid(ncid,'kp',idv_kp)
!     istat = nf_get_var_double(ncid,idv_kp,kp)
!     call fminmaxspv(kp,ndata,fmin,fmax,missing)
!     write(6,"('rdimf: kp min,max=',2f8.2)") fmin,fmax

! f107d:
!     istat = nf_inq_varid(ncid,'f107d',idv_f107d)
!     istat = nf_get_var_double(ncid,idv_f107d,f107d)
!     call fminmaxspv(f107d,ndays,fmin,fmax,missing)
!     write(6,"('rdimf: f107d min,max=',2f8.2)") fmin,fmax

! f107a:
!     istat = nf_inq_varid(ncid,'f107a',idv_f107a)
!     istat = nf_get_var_double(ncid,idv_f107a,f107a)
!     call fminmaxspv(f107a,ndays,fmin,fmax,missing)
!     write(6,"('rdimf: f107a min,max=',2f8.2)") fmin,fmax

      write(6,"(72('-'),/)")
      end subroutine rdimf
!-----------------------------------------------------------------------
      subroutine alloc_imf

! real,dimension(:),target,allocatable ::
!   uttime,yfrac,bx,by,bz,swvel,swden,kp,f107d,f107a

      integer :: istat
      allocate(date(ndata),stat=istat)
      allocate(uttime(ndata),stat=istat)
      allocate(yfrac(ndata),stat=istat)
      allocate(bx(ndata),stat=istat)
      allocate(by(ndata),stat=istat)
      allocate(bz(ndata),stat=istat)
      allocate(swvel(ndata),stat=istat)
      allocate(swden(ndata),stat=istat)
      allocate(kp(ndata),stat=istat)

!     allocate(days(ndays),stat=istat)
!     allocate(f107d(ndays),stat=istat)
!     allocate(f107a(ndays),stat=istat)

      end subroutine alloc_imf
!-----------------------------------------------------------------------
      subroutine printmissing(i0,i1,name)
      integer,intent(in) :: i0,i1
      character(len=*),intent(in) :: name
!
      write(6,"(/,'>>> getimf: Encountered missing ',a,' in imf',
     |  ' data file (missing=',e12.4,'):')") trim(name),missing
      write(6,"('mtime=',3i4)") modeltime(1:3)
      write(6,"('imf data file = ',a)") trim(imf_ncfile)
      write(6,"('index1=',i8,' date(index1)=',f12.4,' uttime(index1)=',
     |  f8.3)") i0,date(i0),uttime(i0)
      write(6,"('index2=',i8,' date(index2)=',f12.4,' uttime(index2)=',
     |  f8.3)") i1,date(i1),uttime(i1)
      write(6,"('bx(index1)=',e12.4,' bx(index2)=',e12.4)")bx(i0),bx(i1)
      write(6,"('by(index1)=',e12.4,' by(index2)=',e12.4)")by(i0),by(i1)
      write(6,"('bz(index1)=',e12.4,' bz(index2)=',e12.4)")bz(i0),bz(i1)
      write(6,"('swvel(index1)=',e12.4,' swvel(index2)=',e12.4)")
     |  swvel(i0),swvel(i1)
      write(6,"('swden(index1)=',e12.4,' swden(index2)=',e12.4)")
     |  swden(i0),swden(i1)

      end subroutine printmissing
!-----------------------------------------------------------------------
      logical function ismissing(x)
      use params_module,only: spval
      implicit none
      real,intent(in) :: x
      ismissing = .false.
      if (x == spval .or. 
     |   (x >= spval-.1*spval .and. x <= spval+.1)) then
        ismissing = .true. 
!       write(6,"('issmissing true: x=',e12.4,' spval=',e12.4)")
!    |    x,spval
      endif
      end function ismissing
!-----------------------------------------------------------------------
      end module imf_module

