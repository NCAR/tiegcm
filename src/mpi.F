!
      module mpi_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Perform message-passing and related operations in distributed memory 
!   system, e.g., AIX.
!
      use params_module
      use fields_module,only: fields_4d,fields_3d,f4d,f3d,nf3d,
     |  nf4d,nf4d_hist,fsechist,dynpot,phim3d,i_poten,foutput,
     |  tlbc,ulbc,vlbc,tlbc_glb,ulbc_glb,vlbc_glb,emphi3d,
     |  emlam3d,emz3d,fzg,tlbc_nm,ulbc_nm,vlbc_nm,tlbc_nm_glb,
     |  ulbc_nm_glb,vlbc_nm_glb
      use hist_module,only: nfsech
      use input_module,only: ntask_lat,ntask_lon
      implicit none
#ifdef MPI
#include <mpif.h>
      integer ::
     |  irstat(MPI_STATUS_SIZE)      ! mpi receive status
#endif
!
! VT means vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
      integer :: 
     |  ntask,      ! number of mpi tasks
     |  ntaski,     ! number of tasks in lon dimension (from input ntask_lat)
     |  ntaskj,     ! number of tasks in lat dimension (from input ntask_lon)
     |  mytid,      ! task id of current mpi task
     |  mytidi,     ! i coord for current task in task table
     |  mytidj,     ! j coord for current task in task table
     |  lat0,lat1,  ! first and last lats for each task
     |  lon0,lon1,  ! first and last lons for each task
     |  mxlon,mxlat ! max of nlons,nlats owned by all tasks
      integer,allocatable :: 
     |  mylats(:),  ! indices of lats I have
     |  mylons(:)   ! indices of lons I have
      integer,allocatable :: 
     |  itask_table(:,:)  ! 2d table of tasks (i,j)
!
! Task type:
      type task
        integer :: mytid     ! task id
        integer :: mytidi    ! task coord in longitude dimension of task table
        integer :: mytidj    ! task coord in latitude  dimension of task table
        integer :: nlats     ! number of latitudes  calculated by this task
        integer :: nlons     ! number of longitudes calculated by this task
        integer :: lat0,lat1 ! first and last latitude  indices
        integer :: lon0,lon1 ! first and last longitude indices
        integer :: mylats(nlat  ) ! lat indices (only lat0->lat1 are meaningful)
        integer :: mylons(nlonp4) ! lon indices (only lon0->lon1 are meaningful)
      end type task
!
! type(task) :: tasks(ntask) will be made available to all tasks
! (so each task has information about all tasks)
!
      type(task),allocatable :: tasks(:)
!
! Timing stats:
      logical :: do_rtc_mpi=.false.
      real :: tsec
      real :: rtc0_gather2root     , rtcmp_gather2root     =0.,
     |        rtc0_bndlats         , rtcmp_bndlats         =0.,
     |        rtc0_bndlats_f2d     , rtcmp_bndlats_f2d     =0.,
     |        rtc0_bndlons         , rtcmp_bndlons         =0.,
     |        rtc0_bndlons_f3d     , rtcmp_bndlons_f3d     =0.,
     |        rtc0_bndlons_f2d     , rtcmp_bndlons_f2d     =0.,
     |        rtc0_polelat         , rtcmp_polelat         =0.,
     |        rtc0_gatherlons_f3d  , rtcmp_gatherlons_f3d  =0.,
     |        rtc0_scatterlons_f3d , rtcmp_scatterlons_f3d =0.,
     |        rtc0_periodic_f4d    , rtcmp_periodic_f4d    =0.,
     |        rtc0_periodic_f3d    , rtcmp_periodic_f3d    =0.,
     |        rtc0_periodic_f2d    , rtcmp_periodic_f2d    =0.,
     |        rtc0_dynpot          , rtcmp_dynpot          =0.
!
! If not an mpi job, this module contains no subroutines.
!
#ifdef MPI
      contains
!-----------------------------------------------------------------------
      subroutine mp_init
!
! Initialize mpi and get ntask and mytid (called from tgcm.F before input):
!
! Local:
      integer :: ier,n,j,i
      integer,dimension(9) :: 
     |  valid_pecounts=(/1,4,8,12,16,24,32,48,64/)
!
      ier = 0
#ifdef VT
      call vttraceoff(ier)
#endif
      call mpi_init(ier)
      if (ier /= 0) then
        write(6,"('>>> Error from mpi_init: ier=',i4)") ier
        call shutdown('mpi_init error')
      endif
#ifdef VT
      call vtsetup()      ! set user defined states and activities
      write(6,"('mp_init: starting vampir tracing..')")
      call vttraceon(ier) ! start vampir tracing
#endif
      call mpi_comm_size(MPI_COMM_WORLD,ntask,ier)
      if (ier /= 0) then
        write(6,"('>>> Error from mpi_comm_size: ier=',i4)") ier
        call shutdown('mpi_comm_size error')
      endif
      if (.not.any(ntask == valid_pecounts)) then
        write(6,"(/,'>>> This model may not run properly or will ',
     |    'perform poorly with nask=',i4,/,' Please either change ',
     |    'valid_pecounts in sub mp_init of source file mpi.F,',/,
     |    ' or use one of these processor counts for which the model ',
     |    'has been validated:',/,(8i8))") ntask,valid_pecounts
        call shutdown('invalid ntask')
      endif

      call mpi_comm_rank(MPI_COMM_WORLD,mytid,ier)
      if (ier /= 0) then
        write(6,"('>>> Error from mpi_comm_rank: ier=',i4)") ier
        call shutdown('mpi_comm_rank error')
      endif
!
! Allocate array of tasks (user defined type(task)):
      allocate(tasks(0:ntask-1),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating ',i4,' tasks')") ntask
        call shutdown('mp_init')
      endif
      write(6,"('mp_init: ntask=',i3,' mytid=',i3)") ntask,mytid
      end subroutine mp_init
!-----------------------------------------------------------------------
      subroutine mp_distribute
!
! Set up 2-d data decomposition in lat,lon. Define structure array 
!   tasks(ntasks). This is broadcast (mpi_alltoall) to all tasks,
!   so every task has access to info for all tasks.
!
! Set by mp_init:
!   ntask:  Number of mpi tasks this run (from mp_init)
!   mytid:  id of current task (from mp_init)
!
! Set by input (ntask_lon,lat are defined on input to this routine):
!   ntask_lon = ntaski: Number of tasks in longitude direction
!   ntask_lat = ntaskj: Number of tasks in latitude direction
!
! Set by this routine:
!   itask_table(ntaski,ntaskj): 2d table of tasks
!   mytidi: i-coord of current task in task table
!   mytidj: j-coord of current task in task table
!   lat0,lat1: starting and ending indices in latitude for current task
!   lon0,lon1: starting and ending indices in longitude for current task
!   mylats,mylons: indices to lats and lons for current tasks
!
!   tasks(ntasks): array of task types.
!   (note mytidi,j,lat0,1,lon0,1,mylats,lons above are redundant with
!    tasks(mytid), eg, tasks(mytid)%mytidi==mytidi)
!
! Local:
      integer :: i,j,n,ier,irank,nj,ni,ncells
!
! itasks_send(len_task_type,ntask) will be used to send tasks(:) info
!   to all tasks (directly passing mpi derived data types is reportedly 
!   not stable, or not available until MPI 2.x).
!
      integer,parameter :: len_task_type = nlat+nlonp4+9
      integer :: itask_packed(len_task_type)
      integer,allocatable :: 
     |  itasks_send(:,:), ! send buffer
     |  itasks_recv(:,:)  ! send buffer
!
! Set ntaski,j from input ntask_lon,ntask_lat:
      ntaski = ntask_lon
      ntaskj = ntask_lat
!
! ntaski*ntaskj must == ntask
! (This was already checked in input, but you never know...)
      if (ntaski*ntaskj /= ntask) then 
        write(6,"('>>> mp_distribute: ntaski*ntaskj must == ntask.')")
        write(6,"('    ntaski=',i3,' ntaskj=',i3,' ntask=',i3)")
     |    ntaski,ntaskj,ntask
        call shutdown('ntaski*j')
      endif
!
! Allocate and set 2d table of tasks:
      allocate(itask_table(-1:ntaski,-1:ntaskj),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating itable: ntaski,j=',2i3)") 
     |    ntaski,ntaskj
        call mp_close
        call shutdown('itask_table')
      endif
      itask_table(:,:) = MPI_PROC_NULL
      irank = 0
      do j =  0,ntaskj-1
        do i = 0,ntaski-1
          itask_table(i,j) = irank
          if (mytid == irank) then
            mytidi = i
            mytidj = j
          endif
          irank = irank+1
        enddo
      enddo  
!
! Print table to stdout:
      write(6,"(/,'ntask=',i3,' ntaski=',i2,' ntaskj=',i2,
     |  ' Task Table:')") ntask,ntaski,ntaskj
      do j=-1,ntaskj
        write(6,"('j=',i3,' itask_table(:,j)=',100i3)")
     |      j,itask_table(:,j)
      enddo
!
! Calculate start and end indices in lon,lat dimensions for each task:
      call distribute_1d(1,nlonp4,ntaski,mytidi,lon0,lon1)
      call distribute_1d(1,nlat  ,ntaskj,mytidj,lat0,lat1)
      nj = lat1-lat0+1 ! number of latitudes for this task
      ni = lon1-lon0+1 ! number of longitudes for this task
      ncells = nj*ni   ! total number of grid cells for this task
!
! Report my stats to stdout:
      write(6,"(/,'mytid=',i3,' mytidi,j=',2i3,' lat0,1=',2i3,' (',i2,
     |  ') lon0,1=',2i3,' (',i2,') ncells=',i4)") 
     |  mytid,mytidi,mytidj,lat0,lat1,nj,lon0,lon1,ni,ncells
!
! Set mylats, mylons:
      allocate(mylats(nj),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating mylats(nj=',i2,'): ier=',i3)") 
     |    nj,ier
      endif
      do j=1,nj
        mylats(j) = lat0+j-1
      enddo
!
      allocate(mylons(ni),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating mylons(ni=',i2,'): ier=',i3)") 
     |    ni,ier
      endif
      do i=1,ni
        mylons(i) = lon0+i-1
      enddo
!
! Define all task structures with current task values 
! (redundant for alltoall):
!
      do n=0,ntask-1
        tasks(n)%mytid  = mytid
        tasks(n)%mytidi = mytidi
        tasks(n)%mytidj = mytidj
        tasks(n)%nlats  = nj
        tasks(n)%nlons  = ni
        tasks(n)%lat0   = lat0
        tasks(n)%lat1   = lat1
        tasks(n)%lon0   = lon0
        tasks(n)%lon1   = lon1
        tasks(n)%mylats(:) = ispval ! init
        tasks(n)%mylons(:) = ispval ! init
        do j=lat0,lat1
          tasks(n)%mylats(j) = j 
        enddo
        do i=lon0,lon1
          tasks(n)%mylons(i) = i 
        enddo
      enddo
!
! Pack tasks(mytid) into itasks_send:
      allocate(itasks_send(len_task_type,0:ntask-1),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating itasks_send: len_task_type=',
     |    i3,' ntask=',i3)") len_task_type,ntask
      endif
      allocate(itasks_recv(len_task_type,0:ntask-1),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating itasks_recv: len_task_type=',
     |    i3,' ntask=',i3)") len_task_type,ntask
      endif
      do n=0,ntask-1
        itasks_send(1,n) = tasks(mytid)%mytid
        itasks_send(2,n) = tasks(mytid)%mytidi
        itasks_send(3,n) = tasks(mytid)%mytidj
        itasks_send(4,n) = tasks(mytid)%nlats
        itasks_send(5,n) = tasks(mytid)%nlons
        itasks_send(6,n) = tasks(mytid)%lat0
        itasks_send(7,n) = tasks(mytid)%lat1
        itasks_send(8,n) = tasks(mytid)%lon0
        itasks_send(9,n) = tasks(mytid)%lon1
        itasks_send(10:10+nlat-1,n) = tasks(mytid)%mylats(:)
        i = 10+nlat
        itasks_send(i:i+nlonp4-1,n) = tasks(mytid)%mylons(:)
      enddo
!
! Send itasks_send and receive itasks_recv:
      call mpi_alltoall(itasks_send,len_task_type,MPI_INTEGER,
     |                  itasks_recv,len_task_type,MPI_INTEGER,
     |                  MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mpi_alltoall to send/recv itasks')
!
! Unpack itasks_recv into tasks(n) and report to stdout:
      do n=0,ntask-1
        tasks(n)%mytid  = itasks_recv(1,n)
        tasks(n)%mytidi = itasks_recv(2,n)
        tasks(n)%mytidj = itasks_recv(3,n)
        tasks(n)%nlats  = itasks_recv(4,n)
        tasks(n)%nlons  = itasks_recv(5,n)
        tasks(n)%lat0   = itasks_recv(6,n)
        tasks(n)%lat1   = itasks_recv(7,n)
        tasks(n)%lon0   = itasks_recv(8,n)
        tasks(n)%lon1   = itasks_recv(9,n)
        tasks(n)%mylats(:) = itasks_recv(10:10+nlat-1,n)
        i = 10+nlat
        tasks(n)%mylons(:) = itasks_recv(i:i+nlonp4-1,n)
!
        if (n==mytid) then
          write(6,"(/,'Task ',i3,':')") n
          write(6,"('tasks(',i3,')%mytid =',i3)") n,tasks(n)%mytid
          write(6,"('tasks(',i3,')%mytidi=',i3)") n,tasks(n)%mytidi
          write(6,"('tasks(',i3,')%mytidj=',i3)") n,tasks(n)%mytidj
          write(6,"('tasks(',i3,')%nlats =',i3)") n,tasks(n)%nlats 
          write(6,"('tasks(',i3,')%nlons =',i3)") n,tasks(n)%nlons 
          write(6,"('tasks(',i3,')%lat0  =',i3)") n,tasks(n)%lat0  
          write(6,"('tasks(',i3,')%lat1  =',i3)") n,tasks(n)%lat1  
          write(6,"('tasks(',i3,')%lon0  =',i3)") n,tasks(n)%lon0  
          write(6,"('tasks(',i3,')%lon1  =',i3)") n,tasks(n)%lon1  
          write(6,"('tasks(',i3,')%mylats=',/,(15i4))") 
     |      n,tasks(n)%mylats(tasks(n)%lat0:tasks(n)%lat1)
          write(6,"('tasks(',i3,')%mylons=',/,(15i4))") 
     |      n,tasks(n)%mylons(tasks(n)%lon0:tasks(n)%lon1)
        endif
      enddo
!
! Release locally allocated space:
      deallocate(itasks_send)
      deallocate(itasks_recv)
!
! All tasks must have at least 4 longitudes:
      do n=0,ntask-1
        if (tasks(n)%nlons < 4) then
          write(6,"(/,'>>> mp_distribute: each task must carry ',
     |      'at least 4 longitudes. task=',i2,' nlons=',i2)")
     |      n,tasks(n)%nlons
          call shutdown('nlons per task')
        endif
      enddo
!
! mxlon,mxlat are maximum number of lons,lats owned by all tasks:
      mxlon = -9999
      do n=0,ntask-1
        if (tasks(n)%nlons > mxlon) mxlon = tasks(n)%nlons
      enddo
      mxlat = -9999
      do n=0,ntask-1
        if (tasks(n)%nlats > mxlat) mxlat = tasks(n)%nlats
      enddo
!
! For debug:
!     call mp_close
!     call shutdown('mp_distribute')
!
      end subroutine mp_distribute
!-----------------------------------------------------------------------
      subroutine distribute_1d(n1,n2,nprocs,myrank,istart,iend)
!
! Distribute work across a 1d vector(n1->n2) to nprocs. 
! Return start and end indices for proc myrank.
!
! Args:
      integer,intent(in)  :: n1,n2,nprocs,myrank
      integer,intent(out) :: istart,iend
!
! Local:
      integer :: lenproc,iremain,n
!
      n = n2-n1+1
      lenproc = n/nprocs
      iremain = mod(n,nprocs)
      istart = n1 + myrank*lenproc + min(myrank,iremain)
      iend = istart+lenproc-1 
      if (iremain > myrank) iend = iend+1
!
      end subroutine distribute_1d
!-----------------------------------------------------------------------
      subroutine mp_gather2root(ixt,type)
!
! Master task must collect all data (except its own) into its prognostic
!   fields f4d(:)%data prior to writing a history. Slaves send their
!   (lon0:lon1,nlevp1,lat0:lat1) data to the root task.
! If type == 'prim', pass prognostics f4d(:)%data for primary histories.
! If type == 'sech', pass diagnostic fsech(:)%data for secondary histories.
!
! Args:
      integer,intent(in) :: ixt
!
! type=='prim' if primary history data (f4d(:)), secondary ('sech')
!   data (fsech(:)) otherwise.
!
      character(len=*),intent(in) :: type
!
! Local:
      integer :: i
!
      if (do_rtc_mpi) call timer(rtc0_gather2root,tsec,'begin')
!
! Gather for either primary and secondary histories.
! If primary, also gather f3d field subdomains zg to global fzg.
      if (trim(type)=='prim') then
         call mp_gather2root_prim(ixt)
         do i=1,nf3d
           if (trim(f3d(i)%short_name)=='ZG') then
             call mp_gather2root_f3d(f3d(i),fzg)
           endif
         enddo
      elseif(trim(type)=='sech') then
         call mp_gather2root_sechist(ixt)
      else
        write(6,"('mp_gather2root: unknown type ',a,' should be prim ',
     |    'or sech')")
        call shutdown('mp_gather2root')
      endif
!
! Gather lower boundary condition arrays tlbc,ulbc,vlbc and 
! tlbc_nm,ulbc_nm,vlbc_nm:
      call mp_gather2root_lbc(ixt)
!
      if (do_rtc_mpi) then
        call timer(rtc0_gather2root,tsec,'end')
        rtcmp_gather2root = rtcmp_gather2root+tsec
      endif
!
      end subroutine mp_gather2root
!-----------------------------------------------------------------------      
      subroutine mp_gather2root_prim(ixt)
!
! Master task must collect all data (except its own) into its prognostic 
!   fields f4d(:)%data prior to writing a history. Slaves send their 
!   (lon0:lon1,nlevp1,lat0:lat1) data to the root task.
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: isrc,ireqrecv,ireqsend,ier,n,i,j,k,l,len,
     |  ilon0,ilon1,ilat0,ilat1,nlons,nlats,msgtag
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable,save :: rcvbuf(:,:,:,:),sndbuf(:,:,:,:)
!
! Allocate send and receive buffers:
! (mxlon,mxlat are max number of lons,lats held by all tasks)
!
! Prognostics (primary history data f4d(i)%data(k,i,j,ixt):
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,mxlon,mxlat,nf4d),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_prim: error allocating rcvbuf:',
     |      ' mxlon=',i3,' nlevp1=',i3,' mxlat=',i3,' nf4d=',i3,
     |      ' ier=',i4)") mxlon,nlevp1,mxlat,nf4d,ier
      endif
!
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevp1,mxlon,mxlat,nf4d),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_prim: error allocating sndbuf:',
     |      ' mxlon=',i3,' nlevp1=',i3,' mxlat=',i3,' nf4d=',i3,
     |      ' ier=',i4)") mxlon,nlevp1,mxlat,nf4d,ier
      endif
      rcvbuf = 0. ; sndbuf = 0.
      len = nlevp1*mxlon*mxlat*nf4d ! buffer length
!
! Root receives from all other tasks:
      if (mytid==0) then ! root task
!
! If primary history, root task defines foutput at its own subdomain:
        do i=1,nf4d
          foutput(:,lon0:lon1,lat0:lat1,i) =
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,ixt)
!         call fminmax(
!    |      foutput(:,lon0:lon1,lat0:lat1,i),
!    |      (lon1-lon0+1)*nlevp1*(lat1-lat0+1),fmin,fmax)
!         write(6,"('Root copied its own domain into foutput: i=',
!    |      i3,' min,max=',2e12.4)") i,fmin,fmax
        enddo
!
! Receive subdomains from slave tasks:
        do n=1,ntask-1   ! receive from slaves only
          msgtag = n
          call mpi_irecv(rcvbuf,len,MPI_REAL8,n,msgtag,MPI_COMM_WORLD,
     |      ireqrecv,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gather2root irecv from n')
!
! Wait for receives to complete:
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) then
            call handle_mpi_err(ier,'mp_gather2root wait for n')
          else
!           write(6,"('mp_gather2root: root received from task ',i3)") n
          endif
!
! Transfer from receive buffer to fields:
          ilon0 = tasks(n)%lon0
          ilon1 = tasks(n)%lon1
          ilat0 = tasks(n)%lat0
          ilat1 = tasks(n)%lat1
          nlons = ilon1-ilon0+1 
          nlats = ilat1-ilat0+1 
          do i=1,nf4d
!
! Root task stores data in full 3d grid array foutput. This is a pointer
! declared in fields module (fields.F), and allocated by the root task 
! only in sub allocdata: foutput(nlevp1,nlonp4,nlat,nf4d).
!
            foutput(:,ilon0:ilon1,ilat0:ilat1,i) =
     |        rcvbuf(:,1:nlons,1:nlats,i)
!           call fminmax(
!    |        foutput(:,ilon0:ilon1,ilat0:ilat1,i),nlons*nlevp1*nlats,
!    |   	fmin,fmax)
!           if (i==1) write(6,"(' ')")
!           write(6,"('mp_gather2root_prim: received field ',i3,
!    |        ' from task ',i3,' fmin,max=',2e12.4)") i,n,fmin,fmax
          enddo ! i=1,nf4d
!
        enddo ! n=1,ntask-1
!
! Non-root tasks send to master:
      else
!
! Load send buffer:
        sndbuf(:,:,:,:) = 0.
        nlons = lon1-lon0+1 
        nlats = lat1-lat0+1 
!
! Prognostic fields:
        do i=1,nf4d
          sndbuf(:,1:nlons,1:nlats,i) = 
     |      f4d(i)%data(:,lon0:lon1,lat0:lat1,ixt)
!         call fminmax(sndbuf(:,1:nlons,1:nlats,i),nlons*nlevp1*nlats,
!    |      fmin,fmax)
!         write(6,"('mp_gather2root_prim: task ',i3,' sending prog ',
!    |      'field ',i2,' to root: fmin,max=',2e12.4)") 
!    |      mytid,i,fmin,fmax
        enddo
!
! Send to root:
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gather2root_prim isend to root')
!
! Wait for send to complete: 
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0) then
          call handle_mpi_err(ier,
     |    'mp_gather2root_prim wait for send to root')
        endif
      endif ! root or slave
      end subroutine mp_gather2root_prim
!-----------------------------------------------------------------------      
      subroutine mp_gather2root_f3d(f3d,fout)
!
! Master task must collect all data (except its own) into its prognostic 
!   fields f4d(:)%data prior to writing a history. Slaves send their 
!   (lon0:lon1,nlevp1,lat0:lat1) data to the root task.
!
! Args:
      type(fields_3d),intent(in) :: f3d
      real,intent(out) :: fout(nlevp1,nlonp4,nlat) ! see allocdata.F
!
! Local:
      integer :: isrc,ireqrecv,ireqsend,ier,n,i,j,k,l,len,
     |  ilon0,ilon1,ilat0,ilat1,nlons,nlats,msgtag
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable,save :: rcvbuf(:,:,:),sndbuf(:,:,:)
!
! Allocate send and receive buffers:
! (mxlon,mxlat are max number of lons,lats held by all tasks)
!
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,mxlon,mxlat),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_f3d: error allocating rcvbuf:',
     |      ' mxlon=',i3,' nlevp1=',i3,' mxlat=',i3,' ier=',i4)") 
     |      mxlon,nlevp1,mxlat,ier
      endif
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevp1,mxlon,mxlat),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_f3d: error allocating sndbuf:',
     |      ' mxlon=',i3,' nlevp1=',i3,' mxlat=',i3,' ier=',i4)") 
     |      mxlon,nlevp1,mxlat,ier
      endif
      rcvbuf = 0. ; sndbuf = 0.
      len = nlevp1*mxlon*mxlat ! buffer length
!
! Root receives from all other tasks:
      if (mytid==0) then ! root task
!
! If primary history, root task defines foutput at its own subdomain:
        fout(:,lon0:lon1,lat0:lat1) =
     |    f3d%data(:,lon0:lon1,lat0:lat1)
!       call fminmax(
!    |    fout(:,lon0:lon1,lat0:lat1),
!    |    (lon1-lon0+1)*nlevp1*(lat1-lat0+1),fmin,fmax)
!       write(6,"('Root copied its own domain into foutput: ',
!    |    ' min,max=',2e12.4)") fmin,fmax
!
! Receive subdomains from slave tasks:
        do n=1,ntask-1   ! receive from slaves only
          msgtag = n
          call mpi_irecv(rcvbuf,len,MPI_REAL8,n,msgtag,MPI_COMM_WORLD,
     |      ireqrecv,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gather2root_f3d irecv from n')
!
! Wait for receives to complete:
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) then
            call handle_mpi_err(ier,'mp_gather2root_f3d wait for n')
          else
!           write(6,"('mp_gather2root_f3d: root received from task ',
!    |        i3)") n
          endif
!
! Transfer from receive buffer to fields:
          ilon0 = tasks(n)%lon0
          ilon1 = tasks(n)%lon1
          ilat0 = tasks(n)%lat0
          ilat1 = tasks(n)%lat1
          nlons = ilon1-ilon0+1 
          nlats = ilat1-ilat0+1 
!
! Root task stores data in full 3d grid array foutput. This is a pointer
! declared in fields module (fields.F), and allocated by the root task 
! only in sub allocdata: foutput(nlevp1,nlonp4,nlat,nf4d).
!
          fout(:,ilon0:ilon1,ilat0:ilat1) = rcvbuf(:,1:nlons,1:nlats)
 
!         call fminmax(
!    |      fout(:,ilon0:ilon1,ilat0:ilat1),nlons*nlevp1*nlats,
!    |        fmin,fmax)
!         write(6,"('mp_gather2root_f3d: received field',
!    |      ' from task ',i3,' fmin,max=',2e12.4)") n,fmin,fmax
 
        enddo ! n=1,ntask-1
!
! Non-root tasks send to master:
      else
!
! Load send buffer:
        sndbuf(:,:,:) = 0.
        nlons = lon1-lon0+1 
        nlats = lat1-lat0+1 
!
! Prognostic fields:
        sndbuf(:,1:nlons,1:nlats) = f3d%data(:,lon0:lon1,lat0:lat1)

!       call fminmax(sndbuf(:,1:nlons,1:nlats),nlons*nlevp1*nlats,
!    |    fmin,fmax)
!       write(6,"('mp_gather2root_f3d: task ',i3,' sending prog ',
!    |    'field to root: fmin,max=',2e12.4)") mytid,fmin,fmax
 
! Send to root:
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gather2root_f3d isend to root')
!
! Wait for send to complete: 
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0) then
          call handle_mpi_err(ier,
     |    'mp_gather2root_f3d wait for send to root')
        endif
      endif ! root or slave
      end subroutine mp_gather2root_f3d
!-----------------------------------------------------------------------
      subroutine mp_gather2root_sechist(ixt)
!
! Master task must collect all data (except its own) into its secondary 
!   fields array fsechist(:)%data prior to writing a history. Slaves send 
!   their data to the root task.
!
! 11/11/05 btf: rewritten for new addfld (this replaces old mp_gather2root_sech)
! 5/27/08 btf: Corrected use of sndbuf and rcvbuf to use 1:nlons,1:nlats
!              instead of ilon0:ilon1,ilat0:ilat1. This was causing slave
!              tasks to hang on deallocation of sndbuf.
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: isrc,ireqrecv,ireqsend,ier,n,i,j,k,l,len,ilon0,ilon1,
     |  ilat0,ilat1,msgtag,mxlev,nlons,nlats
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable,save :: rcvbuf(:,:,:,:),sndbuf(:,:,:,:)
!      
! Allocate send and receive buffers:
! Secondary history data fsechhist(i)%data(lon,lat,lev) (mag or geo)
!
      mxlev = nlevp1
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxlon,mxlat,mxlev,nfsech),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_sech: error allocating rcvbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' mxlev=',i3,' nfsech=',
     |      i3,' ier=',i4)") mxlon,mxlat,mxlev,nfsech,ier
      endif
!
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxlon,mxlat,mxlev,nfsech),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_sech: error allocating sndbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' mxlev=',i3,' nfsech=',
     |      i3,' ier=',i4)") mxlon,mxlat,mxlev,nfsech,ier
      endif
      rcvbuf = 0. ; sndbuf = 0.
      len = mxlon*mxlat*mxlev*nfsech ! buffer length
!
! Root receives from all other tasks:
      if (mytid==0) then !-------------- root task (receive) -------------
!
! Receive subdomains from slave tasks:
! don't process the mag fields, since they are only known by the
! master task - so the gather is not done
!
        do n=1,ntask-1   ! receive from slaves only
          msgtag = n
          call mpi_irecv(rcvbuf,len,MPI_REAL8,n,msgtag,MPI_COMM_WORLD,
     |      ireqrecv,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gather2root_sech irecv from n')
!
! Wait for receives to complete:
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) then
            call handle_mpi_err(ier,'mp_gather2root_sec wait for n')
          else
!           write(6,"('mp_gather2root: root received from task ',i3)") n
          endif
!
! Transfer from receive buffer to fields:
          ilon0 = tasks(n)%lon0
          ilon1 = tasks(n)%lon1
          ilat0 = tasks(n)%lat0
          ilat1 = tasks(n)%lat1
          nlons = ilon1-ilon0+1
          nlats = ilat1-ilat0+1
!
! Load secondary history fields from receive buffer:
          do i = 1,nfsech
            if (.not.fsechist(i)%prognostic.and. 
     |          .not.fsechist(i)%task0_only.and.
     |          associated(fsechist(i)%data)) then
!
! 3d field (lon,lat,lev):
              if (fsechist(i)%ndims==3) then 
                fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,:)=
     |      	          rcvbuf(1:nlons,1:nlats,:,i)
              else 
!
! 2d field (lon,lat):
                fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,1)=
     |      	          rcvbuf(1:nlons,1:nlats,1,i)
              endif ! 2d or 3d
            endif ! not prognostic and not task0 only
          enddo ! i=1,nfsech
        enddo  ! 1,ntask-1 receive loop
!
! Non-root tasks send to master:
      else              !-------------- non-root task (send) -----------
!
! Load send buffer:
        sndbuf(:,:,:,:) = 0.
!
! Load 3d secondary history field subdomains to send buffer:
        nlons = lon1-lon0+1
        nlats = lat1-lat0+1
        do i = 1,nfsech  
          if (.not.fsechist(i)%prognostic.and. 
     |        .not.fsechist(i)%task0_only.and.
     |        associated(fsechist(i)%data)) then
!
! 3d send buffer:
            if (fsechist(i)%ndims==3) then     ! 3d send buffer
                        sndbuf(1:nlons,1:nlats,:,i) =
     |        fsechist(i)%data(lon0:lon1,lat0:lat1,:)
            else                               ! 2d send buffer
!
! 2d send buffer:
                        sndbuf(1:nlons,1:nlats,1,i) =
     |        fsechist(i)%data(lon0:lon1,lat0:lat1,1)
            endif
          endif
        enddo ! i = 1,nfsech
!
! Send to root:
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gather2root isend to root')
!
! Wait for send to complete: 
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0) call handle_mpi_err(ier,
     |    'mp_gather2root wait for send to root')
      endif ! root or slave
      end subroutine mp_gather2root_sechist
!-----------------------------------------------------------------------
      subroutine mp_gather2root_lbc(ixt)
!
! Gather t,u,v lbc at t and t-1 to root task (tlbc,ulbc,vlbc subdomains are gathered
! to tlbc_glb,ulbc_glb,vlbc_glb and tlbc_nm_glb,ulbc_nm_glb,vlbc_nm_glb , see fields.F)
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: ier,len,msgtag,ireqrecv,nlats,nlons,ireqsend,
     |  ilat1,ilon0,n,ilat0,ilon1
      integer :: idest = 0
      real :: fmin,fmax
      real,allocatable,save :: rcvbuf(:,:,:),sndbuf(:,:,:)
!
! Allocate send and receive buffers:
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxlon,mxlat,6),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_lbc: error allocating rcvbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' ier=',i4)") mxlon,mxlat,ier
      endif
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxlon,mxlat,6),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_lbc: error allocating sndbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' ier=',i4)") mxlon,mxlat,ier
      endif
      rcvbuf = 0. ; sndbuf = 0.
      len = mxlon*mxlat*6 ! buffer length
!
! Root receives from all other tasks:
      if (mytid==0) then ! root task
!
! Root task defines output at its own subdomain:
        tlbc_glb(lon0:lon1,lat0:lat1) = tlbc(lon0:lon1,lat0:lat1)
        ulbc_glb(lon0:lon1,lat0:lat1) = ulbc(lon0:lon1,lat0:lat1)
        vlbc_glb(lon0:lon1,lat0:lat1) = vlbc(lon0:lon1,lat0:lat1)
        tlbc_nm_glb(lon0:lon1,lat0:lat1) = tlbc_nm(lon0:lon1,lat0:lat1)
        ulbc_nm_glb(lon0:lon1,lat0:lat1) = ulbc_nm(lon0:lon1,lat0:lat1)
        vlbc_nm_glb(lon0:lon1,lat0:lat1) = vlbc_nm(lon0:lon1,lat0:lat1)
!
! Receive subdomains from slave tasks:
        do n=1,ntask-1   ! receive from slaves only
          msgtag = n
          call mpi_irecv(rcvbuf,len,MPI_REAL8,n,msgtag,MPI_COMM_WORLD,
     |      ireqrecv,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gather2root_lbc irecv from n')
!
! Wait for receives to complete:
          call mpi_wait(ireqrecv,irstat,ier)
          if (ier /= 0) then
            call handle_mpi_err(ier,'mp_gather2root_lbc wait for n')
          else
!           write(6,"('mp_gather2root_lbc: root received from task ',
!    |        i3)") n
          endif
!
! Transfer from receive buffer to fields:
          ilon0 = tasks(n)%lon0
          ilon1 = tasks(n)%lon1
          ilat0 = tasks(n)%lat0
          ilat1 = tasks(n)%lat1
          nlons = ilon1-ilon0+1 
          nlats = ilat1-ilat0+1 
!
! Root task stores data in full 3d grid arrays tlbc_glb,ulbc_glb,vlbc_glb. 
! (see fields.F)
!
          tlbc_glb(ilon0:ilon1,ilat0:ilat1) = rcvbuf(1:nlons,1:nlats,1)
          ulbc_glb(ilon0:ilon1,ilat0:ilat1) = rcvbuf(1:nlons,1:nlats,2)
          vlbc_glb(ilon0:ilon1,ilat0:ilat1) = rcvbuf(1:nlons,1:nlats,3)
          tlbc_nm_glb(ilon0:ilon1,ilat0:ilat1)=rcvbuf(1:nlons,1:nlats,4)
          ulbc_nm_glb(ilon0:ilon1,ilat0:ilat1)=rcvbuf(1:nlons,1:nlats,5)
          vlbc_nm_glb(ilon0:ilon1,ilat0:ilat1)=rcvbuf(1:nlons,1:nlats,6)
        enddo ! n=1,ntask-1
!
! Non-root tasks send to master:
      else
!
! Load send buffer:
        sndbuf(:,:,:) = 0.
        nlons = lon1-lon0+1 
        nlats = lat1-lat0+1 
!
! Define send buffer:
        sndbuf(1:nlons,1:nlats,1) = tlbc(lon0:lon1,lat0:lat1)
        sndbuf(1:nlons,1:nlats,2) = ulbc(lon0:lon1,lat0:lat1)
        sndbuf(1:nlons,1:nlats,3) = vlbc(lon0:lon1,lat0:lat1)
        sndbuf(1:nlons,1:nlats,4) = tlbc_nm(lon0:lon1,lat0:lat1)
        sndbuf(1:nlons,1:nlats,5) = ulbc_nm(lon0:lon1,lat0:lat1)
        sndbuf(1:nlons,1:nlats,6) = vlbc_nm(lon0:lon1,lat0:lat1)
!
! Send to root:
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gather2root_lbc isend to root')
!
! Wait for send to complete: 
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0) then
          call handle_mpi_err(ier,
     |    'mp_gather2root_lbc wait for send to root')
        endif
      endif ! root or slave
      end subroutine mp_gather2root_lbc
!-----------------------------------------------------------------------
      subroutine mp_bndlats(f,mxf,ixt)
!
! Exchange boundary latitude data between tasks. 
! Each task sends its lat0,lat0+1 data to task jprev, and its lat1-1,lat1
!   data to jnext.
! Each task receives its lat0-2,lat0-1 data from task jprev, and its
!   lat1+1,lat1+2 data from task jnext.
!
! Args:
      integer,intent(in) :: ixt,mxf
      type(fields_4d),intent(inout) :: f(mxf) ! assume 4d data in f(n)%data
!
! Local:
      integer :: n,nn,nlons,ier,len,jnext,jprev,nflds,jsend0,jsend1,
     |  jrecv0,jrecv1,lendat
      real,allocatable,save :: 
     |  sndbuf0(:,:,:,:), ! send buffer for lat0  ,lat0+1 (k,i,2,nf)
     |  sndbuf1(:,:,:,:), ! send buffer for lat1-1,lat1   (k,i,2,nf)
     |  rcvbuf0(:,:,:,:), ! recv buffer for lat0-2,lat0-1 (k,i,2,nf)
     |  rcvbuf1(:,:,:,:)  ! recv buffer for lat1+1,lat1+2 (k,i,2,nf)
      real :: fmin,fmax
!
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats','Communication',ier)
      call vtbegin(100,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_bndlats,tsec,'begin')
!
! Allocate send and receive buffers (only fields with %mpi==true):
      nlons = lon1-lon0+1
      nflds = 0
      do n=1,mxf
        if (f(n)%mpi) nflds = nflds+1 
      enddo
      if (nflds==0) then
        write(6,"('>>> WARNING mp_bndlats: no f(:)%mpi are true --',
     |    ' returning.')")
        return
      endif
      if (.not.allocated(sndbuf0)) then
        allocate(sndbuf0(nlevp1,mxlon,2,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlats: error allocating sndbuf0.')")
      endif
      if (.not.allocated(sndbuf1)) then
        allocate(sndbuf1(nlevp1,mxlon,2,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlats: error allocating sndbuf1.')")
      endif
      if (.not.allocated(rcvbuf0)) then
        allocate(rcvbuf0(nlevp1,mxlon,2,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlats: error allocating rcvbuf0.')")
      endif
      if (.not.allocated(rcvbuf1)) then
        allocate(rcvbuf1(nlevp1,mxlon,2,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlats: error allocating rcvbuf1.')")
      endif
      rcvbuf0 = 0. ; sndbuf0 = 0.
      rcvbuf1 = 0. ; sndbuf1 = 0.
      len    = nlevp1*mxlon*2*nflds 
      lendat = nlevp1*mxlon*2
!
! Locate adjacent tasks (includes null tasks top and bottom)
      jprev = itask_table(mytidi,mytidj-1) ! task above (south)
      jnext = itask_table(mytidi,mytidj+1) ! task below (north)
!
! Load sndbuf0 with lat0,lat0+1 and sndbuf1 with lat1-1,lat1:
      nn = 0
      do n=1,mxf
        if (f(n)%mpi) then
          nn = nn+1
          sndbuf0(:,1:nlons,:,nn) = f(n)%data(:,lon0:lon1,lat0:lat0+1,
     |      ixt)
          sndbuf1(:,1:nlons,:,nn) = f(n)%data(:,lon0:lon1,lat1-1:lat1,
     |      ixt)
        endif
      enddo
!
! Send lat0:lat0+1 (sndbuf0) to jprev:
      call mpi_isend(sndbuf0,len,MPI_REAL8,jprev,1,MPI_COMM_WORLD,
     |  jsend0,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats send0 to jprev')
!
! Send lat1-1:lat1 (sndbuf1) to jnext:
      call mpi_isend(sndbuf1,len,MPI_REAL8,jnext,1,MPI_COMM_WORLD,
     |  jsend1,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats send1 to jnext')
!
! Receive lat0-2:lat0-1 (rcvbuf0) from jprev:
      call mpi_irecv(rcvbuf0,len,MPI_REAL8,jprev,1,MPI_COMM_WORLD,
     |  jrecv0,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats recv0 fm jprev')
!
! Receive lat1+1:lat1+2 (rcvbuf1) from jnext:
      call mpi_irecv(rcvbuf1,len,MPI_REAL8,jnext,1,MPI_COMM_WORLD,
     |  jrecv1,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats recv1 fm jnext')
!
! Wait for completions:
      call mpi_wait(jsend0,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats wait for send0')
      call mpi_wait(jsend1,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats wait for send1')
      call mpi_wait(jrecv0,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats wait for recv0')
      call mpi_wait(jrecv1,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats wait for recv1')
!
! Copy lat0-2:lat0-1 from rcvbuf0, and lat1+1:lat1+2 from rcvbuf1:
      nn = 0
      do n=1,mxf
        if (f(n)%mpi) then
          nn = nn+1
          if (lat0 /= 1) then
            f(n)%data(:,lon0:lon1,lat0-2:lat0-1,ixt) = 
     |        rcvbuf0(:,1:nlons,:,nn)
          endif
          if (lat1 /= nlat) then
            f(n)%data(:,lon0:lon1,lat1+1:lat1+2,ixt) = 
     |        rcvbuf1(:,1:nlons,:,nn)
          endif
        endif
      enddo
      if (do_rtc_mpi) then
        call timer(rtc0_bndlats,tsec,'end')
        rtcmp_bndlats = rtcmp_bndlats+tsec
      endif
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats','Communication',ier)
      call vtend(100,ier)
#endif
      end subroutine mp_bndlats
!-----------------------------------------------------------------------
      subroutine mp_bndlats_f2d(f,i0,i1,j0,j1,nf)
!
! Exchange boundary latitude data between tasks. 
! Each task sends its lat0,lat0+1 data to task jprev, and its lat1-1,lat1
!   data to jnext.
! Each task receives its lat0-2,lat0-1 data from task jprev, and its
!   lat1+1,lat1+2 data from task jnext.
!
! Args:
      integer,intent(in) :: i0,i1,j0,j1,nf
      real,intent(inout) :: f(i0-2:i1+2,j0-2:j1+2,nf)
!
! Local:
      integer :: n,nlons,ier,len,jnext,jprev,nflds,jsend0,jsend1,
     |  jrecv0,jrecv1
      real,allocatable,save :: 
     |  sndbuf0(:,:,:), ! send buffer for j0  ,j0+1 (i,2,nf)
     |  sndbuf1(:,:,:), ! send buffer for j1-1,j1   (i,2,nf)
     |  rcvbuf0(:,:,:), ! recv buffer for j0-2,j0-1 (i,2,nf)
     |  rcvbuf1(:,:,:)  ! recv buffer for j1+1,j1+2 (i,2,nf)
      real :: fmin,fmax
      integer :: ireq(4),istat(MPI_STATUS_SIZE,4),nfsave=0
!
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats_f2d','Communication',ier)
      call vtbegin(100,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_bndlats_f2d,tsec,'begin')
!
! Allocate send and receive buffers (only fields with %mpi==true):
! This may be called w/ different nf -- when nf is greater than the
! previous nf, deallocate the buffers so they will be reallocated w/
! the new larger nf (note 1:nf in the isend/irecv calls):
!
      if (nf > nfsave) then
        if (allocated(sndbuf0)) deallocate(sndbuf0)
        if (allocated(sndbuf1)) deallocate(sndbuf1)
        if (allocated(rcvbuf0)) deallocate(rcvbuf0)
        if (allocated(rcvbuf1)) deallocate(rcvbuf1)
        nfsave = nf
      endif
      if (.not.allocated(sndbuf0)) then
        allocate(sndbuf0(mxlon,2,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlats_f2d: error allocating sndbuf0.')")
      endif
      if (.not.allocated(sndbuf1)) then
        allocate(sndbuf1(mxlon,2,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlats_f2d: error allocating sndbuf1.')")
      endif
      if (.not.allocated(rcvbuf0)) then
        allocate(rcvbuf0(mxlon,2,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlats_f2d: error allocating rcvbuf0.')")
      endif
      if (.not.allocated(rcvbuf1)) then
        allocate(rcvbuf1(mxlon,2,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlats_f2d: error allocating rcvbuf1.')")
      endif
      rcvbuf0 = 0. ; sndbuf0 = 0.
      rcvbuf1 = 0. ; sndbuf1 = 0.
      len    = mxlon*2*nf
      nlons = i1-i0+1
!
! Locate adjacent tasks (includes null tasks top and bottom)
      jprev = itask_table(mytidi,mytidj-1) ! task above (south)
      jnext = itask_table(mytidi,mytidj+1) ! task below (north)
!
! Load sndbuf0 with j0,j0+1 and sndbuf1 with j1-1,j1:
      do n=1,nf
        sndbuf0(1:nlons,:,n) = f(i0:i1,j0  :j0+1,n)
        sndbuf1(1:nlons,:,n) = f(i0:i1,j1-1:j1  ,n)
      enddo
!
! Send j0:j0+1 (sndbuf0) to jprev:
      call mpi_isend(sndbuf0(:,:,1:nf),len,MPI_REAL8,jprev,1,
     |  MPI_COMM_WORLD,jsend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_f2d send0 to jprev')
!
! Send j1-1:j1 (sndbuf1) to jnext:
      call mpi_isend(sndbuf1(:,:,1:nf),len,MPI_REAL8,jnext,1,
     |  MPI_COMM_WORLD,jsend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_f2d send1 to jnext')
!
! Receive j0-2:j0-1 (rcvbuf0) from jprev:
      call mpi_irecv(rcvbuf0(:,:,1:nf),len,MPI_REAL8,jprev,1,
     |  MPI_COMM_WORLD,jrecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_f2d recv0 fm jprev')
!
! Receive j1+1:j1+2 (rcvbuf1) from jnext:
      call mpi_irecv(rcvbuf1(:,:,1:nf),len,MPI_REAL8,jnext,1,
     |  MPI_COMM_WORLD,jrecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_f2d recv1 fm jnext')
!
! Wait for completions:
! integer :: ireq(4),istat(4)
      ireq = (/jsend0,jsend1,jrecv0,jrecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlats_f2d waitall')
!
! Copy j0-2:j0-1 from rcvbuf0, and j1+1:j1+2 from rcvbuf1:
      do n=1,nf
        if (j0 /= 1) then
          f(i0:i1,j0-2:j0-1,n) = rcvbuf0(1:nlons,:,n)
        endif
        if (j1 /= nlat) then
          f(i0:i1,j1+1:j1+2,n) = rcvbuf1(1:nlons,:,n)
        endif
      enddo
      if (do_rtc_mpi) then
        call timer(rtc0_bndlats_f2d,tsec,'end')
        rtcmp_bndlats_f2d = rtcmp_bndlats_f2d+tsec
      endif
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats_f2d','Communication',ier)
      call vtend(100,ier)
#endif
      end subroutine mp_bndlats_f2d
!-----------------------------------------------------------------------
      subroutine mp_bndlons(f,mxf,ixt)
!
! Exchange boundary longitude data between tasks.
! Each task sends its lon0,lon0+1 data to task iprev, and its lon1-1,lon1
!   data to inext.
! Each task receives its lon0-2,lon0-1 data from task iprev, and its
!   lon1+1,lon1+2 data from task inext.
! Because of periodic points, lons -1,0 and nlonp4+1,nlonp4+2 are not
!   necessary. Therefore, tasks with mytidi==0 send/recv with dummy iprev 
!   tasks, and tasks with mytidi==ntaski-1 send/recv with dummy inext tasks. 
!   (dummy tasks have mytid==MPI_PROC_NULL, so comm does not take place).
! This routine is called after mp_bndlats, so longitude data passed includes
!   boundary latitudes.
!
! Args:
      integer,intent(in) :: ixt,mxf
      type(fields_4d),intent(inout) :: f(mxf) ! assume 4d data in f(n)%data
!
! Local:
      integer :: n,nn,nlats,ier,len,inext,iprev,nflds,isend0,isend1,
     |  irecv0,irecv1,lendat
      real,allocatable,save ::
     |  sndbuf0(:,:,:,:), ! send buffer for lon0  ,lon0+1 (k,2,j,nf)
     |  sndbuf1(:,:,:,:), ! send buffer for lon1-1,lon1   (k,2,j,nf)
     |  rcvbuf0(:,:,:,:), ! recv buffer for lon0-2,lon0-1 (k,2,j,nf)
     |  rcvbuf1(:,:,:,:)  ! recv buffer for lon1+1,lon1+2 (k,2,j,nf)
      real :: fmin,fmax
#ifdef VT
!     call vtsymdef(101, 'mp_bndlons','Communication',ier)
      call vtbegin(101,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_bndlons,tsec,'begin')
!
! Allocate send and receive buffers (only fields with %mpi==true):
      nlats = (lat1+2)-(lat0-2)+1
      nflds = 0
      do n=1,mxf
        if (f(n)%mpi) nflds = nflds+1 
      enddo
      if (nflds==0) then
        write(6,"('>>> WARNING mp_bndlons: no f(:)%mpi are true --',
     |    ' returning.')")
        return
      endif
      if (.not.allocated(sndbuf0)) then
        allocate(sndbuf0(nlevp1,2,mxlat+4,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlons: error allocating sndbuf0.')")
      endif
      if (.not.allocated(sndbuf1)) then
        allocate(sndbuf1(nlevp1,2,mxlat+4,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlons: error allocating sndbuf1.')")
      endif
      if (.not.allocated(rcvbuf0)) then
        allocate(rcvbuf0(nlevp1,2,mxlat+4,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlons: error allocating rcvbuf0.')")
      endif
      if (.not.allocated(rcvbuf1)) then
        allocate(rcvbuf1(nlevp1,2,mxlat+4,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_bndlons: error allocating rcvbuf1.')")
      endif
      rcvbuf0 = 0. ; sndbuf0 = 0.
      rcvbuf1 = 0. ; sndbuf1 = 0.
      len    = 2*nlevp1*(mxlat+4)*nflds
      lendat = 2*nlevp1*(mxlat+4)
!
! Locate adjacent tasks (includes null tasks on either side):
      iprev = itask_table(mytidi-1,mytidj) ! task to left (west)
      inext = itask_table(mytidi+1,mytidj) ! task to right (east)
!
! Load sndbuf0 with lon0,lon0+1 and sndbuf1 with lon1-1,lon1:
      nn = 0
      do n=1,mxf
        if (f(n)%mpi) then
          nn = nn+1
          sndbuf0(:,:,1:nlats,nn) = 
     |      f(n)%data(:,lon0:lon0+1,lat0-2:lat1+2,ixt)
          sndbuf1(:,:,1:nlats,nn) = 
     |      f(n)%data(:,lon1-1:lon1,lat0-2:lat1+2,ixt)
        endif
      enddo
!
! Send lon0:lon0+1 (sndbuf0) to iprev:
      call mpi_isend(sndbuf0,len,MPI_REAL8,iprev,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons send0 to iprev')
!
! Send lon1-1:lon1 (sndbuf1) to inext:
      call mpi_isend(sndbuf1,len,MPI_REAL8,inext,1,MPI_COMM_WORLD,
     |  isend1,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons send1 to inext')
!
! Receive lon0-2:lon0-1 (rcvbuf0) from iprev:
      call mpi_irecv(rcvbuf0,len,MPI_REAL8,iprev,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons recv0 fm iprev')
!
! Receive lon1+1:lon1+2 (rcvbuf1) from inext:
      call mpi_irecv(rcvbuf1,len,MPI_REAL8,inext,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons recv1 fm inext')
!
! Wait for completions:
      call mpi_wait(isend0,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons wait for send0')
      call mpi_wait(isend1,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons wait for send1')
      call mpi_wait(irecv0,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons wait for recv0')
      call mpi_wait(irecv1,irstat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons wait for recv1')
!
! Copy lon0-2:lon0-1 from rcvbuf0, and lon1+1:lon1+2 from rcvbuf1:
      nn = 0
      do n=1,mxf
        if (f(n)%mpi) then
          nn = nn+1
          if (lon0 /= 1) then
            f(n)%data(:,lon0-2:lon0-1,lat0-2:lat1+2,ixt) = 
     |        rcvbuf0(:,:,1:nlats,nn)
          endif
          if (lon1 /= nlonp4) then
            f(n)%data(:,lon1+1:lon1+2,lat0-2:lat1+2,ixt) = 
     |        rcvbuf1(:,:,1:nlats,nn)
          endif
        endif
      enddo
      if (do_rtc_mpi) then
        call timer(rtc0_bndlons,tsec,'end')
        rtcmp_bndlons = rtcmp_bndlons+tsec
      endif
#ifdef VT
!     call vtsymdef(101, 'mp_bndlons','Communication',ier)
      call vtend(101,ier)
#endif
      end subroutine mp_bndlons
!-----------------------------------------------------------------------
      subroutine mp_bndlons_f3d(f,id1,i0,i1,j0,j1,nf,iprint)
!
!
! Exchange boundary longitude data of 3-d fields between tasks.
! Each task sends its i0,i0+1 data to task iprev, and its i1-1,i1
!   data to inext.
! Each task receives its i0-2,i0-1 data from task iprev, and its
!   i1+1,i1+2 data from task inext.
! Because of periodic points, lons -1,0 and nlonp4+1,nlonp4+2 are not
!   necessary. Therefore, tasks with mytidi==0 send/recv with dummy iprev 
!   tasks, and tasks with mytidi==ntaski-1 send/recv with dummy inext tasks. 
!   (dummy tasks have mytid==MPI_PROC_NULL, so comm does not take place).
! Note that f must be dimensioned (id1,i0-2:i1+2,j0:j1).
!
! Args:
      integer,intent(in) :: id1,i0,i1,j0,j1,nf,iprint
      real,intent(inout) :: f(id1,i0-2:i1+2,j0-2:j1+2,nf)
!
! Local:
      integer :: n,ier,len,inext,iprev,isend0,isend1,nlats,j,
     |  irecv0,irecv1,lendat,jbeg,jend
      integer,save :: nfsave=0, ncalls=0
      real,allocatable,save :: 
     |  sndbuf0(:,:,:,:), ! send buffer for i0  ,i0+1 (id1,2,mxlat+4,nf)
     |  sndbuf1(:,:,:,:), ! send buffer for i1-1,i1   (id1,2,mxlat+4,nf)
     |  rcvbuf0(:,:,:,:), ! recv buffer for i0-2,i0-1 (id1,2,mxlat+4,nf)
     |  rcvbuf1(:,:,:,:)  ! recv buffer for i1+1,i1+2 (id1,2,mxlat+4,nf)
      integer :: ireq(4),istat(MPI_STATUS_SIZE,4)
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f3d','Communication',ier)
      call vtbegin(102,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_bndlons_f3d,tsec,'begin')
!
! Allocate send and receive buffers.
! This may be called w/ different nf -- when nf is greater than the
! previous nf, deallocate the buffers so they will be reallocated w/
! the new larger nf (note 1:nf in the isend/irecv calls):
!
      if (nf > nfsave) then
        if (iprint > 0)
     |  write(6,"('mp_bndlons_f3d: nf=',i3,' nfsave=',i3,' -- will ',
     |    ' reallocate buffers..')") nf,nfsave
        if (allocated(sndbuf0)) deallocate(sndbuf0)
        if (allocated(sndbuf1)) deallocate(sndbuf1)
        if (allocated(rcvbuf0)) deallocate(rcvbuf0)
        if (allocated(rcvbuf1)) deallocate(rcvbuf1)
        nfsave = nf
      endif
      if (.not.allocated(sndbuf0)) then
        allocate(sndbuf0(id1,2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f3d: error allocating sndbuf0.')")
      endif
      if (.not.allocated(sndbuf1)) then
        allocate(sndbuf1(id1,2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f3d: error allocating sndbuf1.')")
      endif
      if (.not.allocated(rcvbuf0)) then
        allocate(rcvbuf0(id1,2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f3d: error allocating rcvbuf0.')")
      endif
      if (.not.allocated(rcvbuf1)) then
        allocate(rcvbuf1(id1,2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f3d: error allocating rcvbuf1.')")
      endif
      rcvbuf0 = 0. ; sndbuf0 = 0.
      rcvbuf1 = 0. ; sndbuf1 = 0.
      nlats = (j1+2)-(j0-2)+1
      len    = id1*2*(mxlat+4)*nf
      lendat = id1*2*(mxlat+4)
      if (iprint > 0)
     |write(6,"('mp_bndlons_f3d: ncalls=',i6,' id1=',i4,' mxlat+4=',i4,
     |  ' nf=',i4,' len=',i6,' lendat=',i6)") ncalls,id1,mxlat+4,nf,len,
     |  lendat
!
! Locate adjacent tasks (includes null tasks on either side):
      iprev = itask_table(mytidi-1,mytidj) ! task to left (west)
      inext = itask_table(mytidi+1,mytidj) ! task to right (east)
      jbeg = j0-2
      jend = j1+2
!
! Load sndbuf0 with i0,i0+1 and sndbuf1 with i1-1,i1:
      do n=1,nf
        do j=jbeg,jend
          sndbuf0(:,:,j-jbeg+1,n) = f(:,i0:i0+1,j,n)
          sndbuf1(:,:,j-jbeg+1,n) = f(:,i1-1:i1,j,n)
          if (iprint > 0)
     |    write(6,"('mp_bndlons_f3d: n=',i3,' j0,1=',2i4,' j=',i3,
     |      ' sndbuf0(id1-1,:,j-j0+1,n)=',2e12.4,
     |      ' sndbuf1(id1-1,:,j-j0+1,n)=',2e12.4)") n,j0,j1,j,
     |      sndbuf0(id1-1,:,j-j0+1,n),sndbuf1(id1-1,:,j-j0+1,n)
        enddo
      enddo
!
! Send i0:i0+1 (sndbuf0) to task iprev:
      call mpi_isend(sndbuf0(:,:,:,1:nf),len,MPI_REAL8,iprev,1,
     |  MPI_COMM_WORLD,isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f3d send0 to iprev')
!
! Send i1-1:i1 (sndbuf1) to task inext:
      call mpi_isend(sndbuf1(:,:,:,1:nf),len,MPI_REAL8,inext,1,
     |  MPI_COMM_WORLD,isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f3d send1 to inext')
!
! Receive i0-2:i0-1 (rcvbuf0) from task iprev:
      call mpi_irecv(rcvbuf0(:,:,:,1:nf),len,MPI_REAL8,iprev,1,
     |  MPI_COMM_WORLD,irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f3d recv0 fm iprev')
!
! Receive i1+1:i1+2 (rcvbuf1) from task inext:
      call mpi_irecv(rcvbuf1(:,:,:,1:nf),len,MPI_REAL8,inext,1,
     |  MPI_COMM_WORLD,irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f3d recv1 fm inext')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f3d mpi_waitall')
!
! Copy i0-2:i0-1 from rcvbuf0, and i1+1:i1+2 from rcvbuf1:
      do n=1,nf
        if (i0 /= 1) then
          do j=jbeg,jend
            f(:,i0-2:i0-1,j,n) = rcvbuf0(:,:,j-jbeg+1,n)
            if (iprint > 0)
     |      write(6,"('mp_bndlons_f3d: n=',i3,' j0,1=',2i4,' j=',i3,
     |        ' i0,i1=',2i5,' f(id1-1,i0-2:i0-1,j,n)=',2e12.4)")
     |        n,j0,j1,j,i0,i1,f(id1-1,i1+1:i0-1,j,n)
          enddo
        endif
        if (i1 /= nlonp4) then
          do j=jbeg,jend
            f(:,i1+1:i1+2,j,n) = rcvbuf1(:,:,j-jbeg+1,n)
            if (iprint > 0)
     |      write(6,"('mp_bndlons_f3d: n=',i3,' j0,1=',2i4,' j=',i3,
     |        ' i0,i1=',2i5,' f(id1-1,i1+1:i1+2,j,n)=',2e12.4)")
     |        n,j0,j1,j,i0,i1,f(id1-1,i1+1:i1+2,j,n)
          enddo
        endif
      enddo
      if (do_rtc_mpi) then
        call timer(rtc0_bndlons_f3d,tsec,'end')
        rtcmp_bndlons_f3d = rtcmp_bndlons_f3d+tsec
      endif
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f3d','Communication',ier)
      call vtend(102,ier)
#endif
      end subroutine mp_bndlons_f3d
!-----------------------------------------------------------------------
      subroutine mp_bndlons_f2d(f,i0,i1,j0,j1,nf)
!
! Exchange boundary longitude data of 2-d fields between tasks.
! Each task sends its i0,i0+1 data to task iprev, and its i1-1,i1
!   data to inext.
! Each task receives its i0-2,i0-1 data from task iprev, and its
!   i1+1,i1+2 data from task inext.
! Because of periodic points, lons -1,0 and nlonp4+1,nlonp4+2 are not
!   necessary. Therefore, tasks with mytidi==0 send/recv with dummy iprev 
!   tasks, and tasks with mytidi==ntaski-1 send/recv with dummy inext tasks. 
!   (dummy tasks have mytid==MPI_PROC_NULL, so comm does not take place).
! Note that f must be dimensioned (i0-2:i1+2,j0:j1).
!
! E.g., as called from duv.F:
!     real,dimension(lon0-2:lon1+2,lat0:lat1) ::
!    |  ulbc ! lat diffs (includes lon boundaries) (s8)
!     call mp_bndlons_f2d(ulbc,lon0,lon1,lat0,lat1,1)
!
! Args:
      integer,intent(in) :: i0,i1,j0,j1,nf
      real,intent(inout) :: f(i0-2:i1+2,j0-2:j1+2,nf)
!
! Local:
      integer :: n,ier,len,inext,iprev,isend0,isend1,nlats,j,
     |  irecv0,irecv1,jbeg,jend,nfsave=0
      real,allocatable,save :: 
     |  sndbuf0(:,:,:), ! send buffer for i0  ,i0+1 (2,mxlat+4,nf)
     |  sndbuf1(:,:,:), ! send buffer for i1-1,i1   (2,mxlat+4,nf)
     |  rcvbuf0(:,:,:), ! recv buffer for i0-2,i0-1 (2,mxlat+4,nf)
     |  rcvbuf1(:,:,:)  ! recv buffer for i1+1,i1+2 (2,mxlat+4,nf)
      real :: fmin,fmax
      integer :: ireq(4),istat(MPI_STATUS_SIZE,4)
      integer :: j0_iprev,j1_iprev,j0_inext,j1_inext
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f2d','Communication',ier)
      call vtbegin(102,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_bndlons_f2d,tsec,'begin')
! 
! Allocate send and receive buffers.
! This may be called w/ different nf -- when nf is greater than the
! previous nf, deallocate the buffers so they will be reallocated w/
! the new larger nf (note 1:nf in the isend/irecv calls):
!
      if (nf > nfsave) then 
        if (allocated(sndbuf0)) deallocate(sndbuf0)
        if (allocated(sndbuf1)) deallocate(sndbuf1)
        if (allocated(rcvbuf0)) deallocate(rcvbuf0)
        if (allocated(rcvbuf1)) deallocate(rcvbuf1)
        nfsave = nf  
      endif
      if (.not.allocated(sndbuf0)) then
        allocate(sndbuf0(2,mxlat+4,nf),stat=ier)
        if (ier /= 0)  
     |    write(6,"('>>> mp_bndlons_f2d: error allocating sndbuf0.')")
      endif
      if (.not.allocated(sndbuf1)) then
        allocate(sndbuf1(2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f2d: error allocating sndbuf1.')")
      endif
      if (.not.allocated(rcvbuf0)) then
        allocate(rcvbuf0(2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f2d: error allocating rcvbuf0.')")
      endif
      if (.not.allocated(rcvbuf1)) then
        allocate(rcvbuf1(2,mxlat+4,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_bndlons_f2d: error allocating rcvbuf1.')")
      endif
      rcvbuf0 = 0. ; sndbuf0 = 0.
      rcvbuf1 = 0. ; sndbuf1 = 0.
      len    = 2*(mxlat+4)*nf
!
! Locate adjacent tasks (includes null tasks on either side):
      iprev = itask_table(mytidi-1,mytidj) ! task to left (west)
      inext = itask_table(mytidi+1,mytidj) ! task to right (east)
      jbeg = j0-2
      jend = j1+2
!
! Load sndbuf0 with i0,i0+1 and sndbuf1 with i1-1,i1:
      do n=1,nf
        do j=jbeg,jend
          sndbuf0(:,j-jbeg+1,n) = f(i0:i0+1,j,n)
          sndbuf1(:,j-jbeg+1,n) = f(i1-1:i1,j,n)
        enddo
      enddo
!
! Send i0:i0+1 (sndbuf0) to task iprev:
      call mpi_isend(sndbuf0(:,:,1:nf),len,MPI_REAL8,iprev,1,
     |  MPI_COMM_WORLD,isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f2d send0 to iprev')
!
! Send i1-1:i1 (sndbuf1) to task inext:
      call mpi_isend(sndbuf1(:,:,1:nf),len,MPI_REAL8,inext,1,
     |  MPI_COMM_WORLD,isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f2d send1 to inext')
!
! Receive i0-2:i0-1 (rcvbuf0) from task iprev:
      call mpi_irecv(rcvbuf0(:,:,1:nf),len,MPI_REAL8,iprev,1,
     |  MPI_COMM_WORLD,irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f2d recv0 fm iprev')
!
! Receive i1+1:i1+2 (rcvbuf1) from task inext:
      call mpi_irecv(rcvbuf1(:,:,1:nf),len,MPI_REAL8,inext,1,
     |  MPI_COMM_WORLD,irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_f2d recv1 fm inext')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier,'mp_bndlons_f2d waitall')
!
! Copy i0-2:i0-1 from rcvbuf0, and i1+1:i1+2 from rcvbuf1:
      do n=1,nf
        if (i0 /= 1) then
          do j=jbeg,jend
            f(i0-2:i0-1,j,n) = rcvbuf0(:,j-jbeg+1,n)
          enddo
        endif
        if (i1 /= nlonp4) then
          do j=jbeg,jend
            f(i1+1:i1+2,j,n) = rcvbuf1(:,j-jbeg+1,n)
          enddo
        endif
      enddo
!
      if (do_rtc_mpi) then
        call timer(rtc0_bndlons_f2d,tsec,'end')
        rtcmp_bndlons_f2d = rtcmp_bndlons_f2d+tsec
      endif
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f3d','Communication',ier)
      call vtend(102,ier)
#endif
      end subroutine mp_bndlons_f2d
!-----------------------------------------------------------------------
      subroutine mp_polelats(ixt)
!
! Make latitudes j=-1,0 from j=2,1, and j=nlat+1,nlat+2 from nlat-1,nlat
! Longitude data are swapped across the poles, with a sign change if
!   necessary (i.e., U or V).
! 3/02: this is an attempt to replace mp_pole_lats with a faster algorithm. 
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: i,ii,j0,j1,if,it,iit,i0,i1,itask,ier,lenrecv,lensend,
     |  nlonhalf,
     |  irecv(ntask), itrecv,
     |  isend(ntask), itsend,
     |  ntrecv,   ! number of tasks from which to receive needlons
     |  ntsend    ! number of tasks to which to send havelons
      integer :: 
     |  newlons(nlonp4),            ! longitudes needed across the pole
     |  lonseast(nlonp4),lonswest(nlonp4),
     |  needlons(nlonp4,0:ntask-1), ! longitudes needed by each task
     |  irecvlons(nlonp4),
     |  isendto(2,0:ntask-1),       ! i0,i1 and tasks to send my lons to
     |  irecvfm(2,0:ntask-1)        ! i0,i1 and tasks to recv my needlons from
      real,allocatable,save ::
     |  rcvbuf(:,:,:,:,:), ! recv for new lats on this task (k,ineed,2,nf,nt)
     |  sndbuf(:,:,:,:,:)  ! send for new lats on remote task (k,ihave,2,nf,nt)
      real :: fmin,fmax,signchange
!
!#ifdef MPI
!
! Serial run of an MPI job (1 task):
      if (ntask==1) then
        call mk_polelat( 0,1,ixt)
        call mk_polelat(-1,2,ixt)
        call mk_polelat(nlat+1,nlat  ,ixt)
        call mk_polelat(nlat+2,nlat-1,ixt)
        return
      endif ! ntask==1
!#else
!
! Serial job (non-MPI):
!     call mk_polelat( 0,1,ixt)
!     call mk_polelat(-1,2,ixt)
!     call mk_polelat(nlat+1,nlat  ,ixt)
!     call mk_polelat(nlat+2,nlat-1,ixt)
!     return
!#endif
!
      if (mytidj /= 0 .and. mytidj /= ntaskj-1) return
!
#ifdef VT
!     call vtsymdef(103, 'mp_polelat','Communication',ier)
      call vtbegin(103,ier)
#endif
!
! Determine longitudes at which east and west hemispheres are exchanged.
! Data at newlons(i) are replaced by data at lonswest(i) for west hemisphere,
!   or lonseast(i) for east hemisphere.
! When nlonp4==76, newlons(:) will be 37-74 and 3-40, i.e.,
!   (1-38) <- (37-74) and (39-76) <- (3-40)
! Another way of looking at it:
!   newlons(1) -> [east,west]lons
!    37-74     ->   1-38 west
!     3-40     ->  39-76 east
!
      nlonhalf = nlonp4/2
      do i=1,nlonhalf
        newlons(i) = nlonhalf+i-2
        lonswest(newlons(i)) = i
      enddo 
      do i=nlonhalf+1,nlonp4
        newlons(i) = 2+i-nlonhalf
        lonseast(newlons(i)) = i
      enddo
!     write(6,"(/,'newlons(:) = ',/,(15i4))") newlons(:)
!     write(6,"('lonswest(:) = ',/,(15i4))") lonswest(:)
!     write(6,"('lonseast(:) = ',/,(15i4))") lonseast(:)
!
! Determine which longitudes each task needs to make new polar boundary 
! latitudes (i.e., across the pole from its own longitudes):
!
      needlons(:,:) = 0
      do it=0,ntask-1
        i0 = tasks(it)%lon0
        i1 = tasks(it)%lon1
        do i=i0,i1
          needlons(i,it) = newlons(i)
        enddo
!       write(6,"('task ',i2,' needlons=',/,(10i4))") it,needlons(:,it)
      enddo
!     write(6,"(/,'My needlons = ',/,(10i4))") needlons(:,mytid) 
!
      irecvfm(1,:) = 9999
      irecvfm(2,:) = -1
      isendto(1,:) = 9999
      isendto(2,:) = -1
      do i=0,ntaski-1 ! loop over tasks in my latitude row in the task table
        itask = tasks(itask_table(i,mytidj))%mytid
        i0 = tasks(itask)%lon0
        i1 = tasks(itask)%lon1
!       if (itask /= mytid) then ! current task is not me
!
! Determine tasks I need to receive my needlons from:
! (1st dimension of irecvfm is 2, for lon0,lon1 to receive)
          do ii=1,nlonp4
            if (needlons(ii,mytid) > 0) then ! I need this lon
              if (needlons(ii,mytid) >= i0 .and. 
     |            needlons(ii,mytid) <= i1) then
                if (needlons(ii,mytid) < irecvfm(1,itask))
     |            irecvfm(1,itask) = needlons(ii,mytid)
                if (needlons(ii,mytid) > irecvfm(2,itask))
     |            irecvfm(2,itask) = needlons(ii,mytid)
              endif
            endif
          enddo
!
! Determine tasks I need to send my longitudes to:
! (1st dimension of isendto is 2, for lon0,lon1 to send)
          do ii=1,nlonp4
            if (needlons(ii,itask) >= lon0 .and. 
     |          needlons(ii,itask) <= lon1) then ! itask needs this lon from me
              if (needlons(ii,itask) < isendto(1,itask))
     |          isendto(1,itask) = needlons(ii,itask)
              if (needlons(ii,itask) > isendto(2,itask))
     |          isendto(2,itask) = needlons(ii,itask)
            endif
          enddo ! ii=1,nlonp4
!       endif ! current task is not me
      enddo ! i=0,ntaski-1  
      where(irecvfm==9999) irecvfm = -1
      where(isendto==9999) isendto = -1
!
! Determine how many tasks will be sent to and received from:
!
      ntrecv = 0 ; ntsend = 0 
      do it=0,ntask-1
        if (irecvfm(1,it) > -1 .and. irecvfm(2,it) > -1) then
          ntrecv = ntrecv+1
        endif
        if (isendto(1,it) > -1 .and. isendto(2,it) > -1) then
          ntsend = ntsend+1  
        endif
      enddo ! it=1,ntask-1
!
! Allocate send buffer.
! sndbuf(:,:,:,:,:), ! send for new lats on remote task (k,i,2,nf,ntsend)
!
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevp1,mxlon,2,nf4d,ntsend),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_polelats: error allocating sndbuf: ',
     |      ' mxlon=',i3,' ntsend=',i3)") mxlon,ntsend
      endif
      sndbuf = 0.
      lensend = nlevp1*mxlon*2*nf4d
!
! Allocate receive buffer.
! rcvbuf(:,:,:,:,:), ! recv for new lats on this task (k,i,2,nf,ntrecv)
!
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,mxlon,2,nf4d,ntrecv),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_polelats: error allocating rcvbuf: ',
     |      ' mxlon=',i3,' ntrecv=',i3)") mxlon,ntrecv
      endif
      rcvbuf = 0.
      lenrecv = nlevp1*mxlon*2*nf4d
!
! Load send buffer:
      j0 = 1 ; j1 = 2              ! Make shem lats 0,-1 from 1,2
      if (mytidj == ntaskj-1) then ! Make nhem lats nlat+1,nlat+2 from nlat-1,nlat
        j0 = nlat-1 ; j1 = nlat
      endif
      iit = 0
      do it=0,ntask-1
        if (isendto(1,it)/=-1.and.isendto(2,it)/=-1) then 
          iit = iit+1
          do if=1,nf4d
            ii = 0
            do i=isendto(1,it),isendto(2,it)
              ii = ii+1 
              sndbuf(:,ii,1:2,if,iit) = f4d(if)%data(:,i,j0:j1,ixt)
            enddo
          enddo ! if=1,nf4d
        endif ! isendto /= -1
      enddo ! it=1,ntaski-1
!
! Post sends of longitudes I have that are needed on other tasks:
!
      itrecv = 0
      itsend = 0
      iit = 0
      do it=0,ntask-1
        if (isendto(1,it)/=-1.and.isendto(2,it)/=-1) then 
          iit = iit+1
          itsend = itsend+1
          call mpi_isend(sndbuf(1,1,1,1,iit),lensend,MPI_REAL8,it,1,
     |      MPI_COMM_WORLD,isend(iit),ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_polelats send')
        endif
      enddo
!
! Post receives of longitudes I need:
!
      iit = 0
      do it=0,ntask-1
        if (irecvfm(1,it)/=-1.and.irecvfm(2,it)/=-1) then 
          iit = iit+1
          call mpi_irecv(rcvbuf(1,1,1,1,iit),lenrecv,MPI_REAL8,it,1,
     |      MPI_COMM_WORLD,irecv(iit),ier)
          if (ier /= 0) call handle_mpi_err(ier,'mp_polelats recv')
        endif
      enddo
!
! Wait for completion of sends:
      do it=1,ntsend
        call mpi_wait(isend(it),irstat,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_polelats wait for send')
      enddo
!
! Wait for completion of receives:
      do it=1,ntrecv
        call mpi_wait(irecv(it),irstat,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_polelats wait for recv')
      enddo
!
! Define f4d at new latitudes from recvbuf:
!
      j0 = -1 ; j1 = 0             ! Receive data for new S lats 0,-1
      if (mytidj == ntaskj-1) then ! Receive data for new N lats nlat+1,nlat+2
        j0 = nlat+1 ; j1 = nlat+2
      endif
      iit = 0
      do it=0,ntask-1
        if (irecvfm(1,it)/=-1.and.irecvfm(2,it)/=-1) then 
          iit = iit+1
          field_loop: do if=1,nf4d
            if (f4d(if)%polesign==0.) cycle field_loop ! n2d,ne,o2p
            ii = 0
            do i=irecvfm(1,it),irecvfm(2,it) ! longitudes that were sent
              ii = ii+1 
!
! Determine if lon to be replaced is east or west:
              if (lon1 <= nlonhalf) then
                irecvlons(ii) = lonswest(i)
              elseif (lon0 > nlonhalf) then
                irecvlons(ii) = lonseast(i)
              else ! task is split between east and west
                irecvlons(ii) = lonswest(i)
                if (i < nlonhalf) irecvlons(ii) = lonseast(i)
              endif
!
! polesign == -1. for u,v velocities (is +1. for all other fields):
              f4d(if)%data(:,irecvlons(ii),j0,ixt)=
     |          rcvbuf(:,ii,2,if,iit)*f4d(if)%polesign
              f4d(if)%data(:,irecvlons(ii),j1,ixt)=
     |          rcvbuf(:,ii,1,if,iit)*f4d(if)%polesign
            enddo
          enddo field_loop
        endif
      enddo
#ifdef VT
!     call vtsymdef(103, 'mp_polelat','Communication',ier)
      call vtend(103,ier)
#endif
      end subroutine mp_polelats
!-----------------------------------------------------------------------
      subroutine mp_gatherlons_f3d(f,k0,k1,i0,i1,j0,j1,nflds)
!
! Gather longitude data in a row of tasks to leftmost task in the row.
! On entry f(k0:k1,i0:i1,j0:j1,nflds) is defined for current task.
! On exit f(k0:k1,nlonp4,j0:j1,nflds) is defined for task with mytidi==0.
!
! Args:
      integer,intent(in) :: k0,k1,i0,i1,j0,j1,nflds
      real,intent(inout) :: f(k0:k1,nlonp4,j0:j1,nflds)
!
! Local:
      integer :: i,k,j,n,nlons,nlats,nlonrecv,nlevs,len,idest,isrc,ier,
     |  isend,irecv,itask,lonrecv0,lonrecv1,mtag
      real :: fmin,fmax
      real,allocatable,save ::
     |  sndbuf(:,:,:,:), ! send buffer (nlevs,mxlon,mxlat,nflds)
     |  rcvbuf(:,:,:,:)  ! recv buffer (nlevs,mxlon,mxlat,nflds)
#ifdef VT
!     call vtsymdef(106, 'mp_gatherlons_f3d','Communication',ier)
      call vtbegin(106,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_gatherlons_f3d,tsec,'begin')
!
! Allocate send and receive buffers:
      nlons = i1-i0+1
      nlevs = k1-k0+1
      nlats = j1-j0+1
      if (nflds==0) then
        write(6,"('>>> WARNING mp_gatherlons_f3d: no fields? nflds=',
     |    i2,' -- returning.')") nflds
        return
      endif
!     if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_gatherlons_f3d: error allocating sndbuf.')")
!     endif
!     if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0) 
     |    write(6,"('>>> mp_gatherlons_f3d: error allocating rcvbuf.')")
!     endif
      sndbuf = 0. ; rcvbuf = 0.
      len = nlevs*mxlon*mxlat*nflds
!
! If mytidi==0, receive from other tasks in my row (mytidi>0,mytidj):
      if (mytidi == 0) then
        do itask=1,ntaski-1
          isrc = itask_table(itask,mytidj)
          mtag = isrc+mytid
          call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,mtag,MPI_COMM_WORLD,
     |      irecv,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gatherlons_f3d recv fm isrc')
          call mpi_wait(irecv,irstat,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_gatherlons_f3d wait for recv0')
!
! Copy data from receive buffer:
          lonrecv0 = tasks(isrc)%lon0
          lonrecv1 = tasks(isrc)%lon1
          nlonrecv = lonrecv1-lonrecv0+1
          do n=1,nflds
            do j=j0,j1
              f(:,lonrecv0:lonrecv1,j,n) = rcvbuf(:,1:nlonrecv,j-j0+1,n)
            enddo ! j=j0,j1
          enddo ! n=1,nflds
        enddo ! itask=1,ntaski-1
!
! If mytidi > 0, load send buffer, and send to task (0,mytidj):
      else ! mytidi /= 0
        idest = itask_table(0,mytidj)
        do n=1,nflds
          do j=j0,j1
            sndbuf(:,1:nlons,j-j0+1,n) = f(:,i0:i1,j,n)
          enddo ! j=j0,j1
        enddo ! n=1,nflds
        mtag = idest+mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,mtag,MPI_COMM_WORLD,
     |    isend,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gatherlons_f3d send0 to idest')
        call mpi_wait(isend,irstat,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_gatherlons_f3d wait for send0')
      endif ! mytidi==0
!
! Free local buffer space:
      deallocate(sndbuf)
      deallocate(rcvbuf)
!
      if (do_rtc_mpi) then
        call timer(rtc0_gatherlons_f3d,tsec,'end')
        rtcmp_gatherlons_f3d = rtcmp_gatherlons_f3d+tsec
      endif
#ifdef VT
!     call vtsymdef(106, 'mp_gatherlons_f3d','Communication',ier)
      call vtend(106,ier)
#endif
      end subroutine mp_gatherlons_f3d
!-----------------------------------------------------------------------
      subroutine mp_scatterlons_f3d(f,k0,k1,i0,i1,j0,j1,nflds,name)
!
! Redistribute longitudes from left most task in j-row to other tasks
!   in the row. Note that nlats is the same for all tasks in each row,
!   but nlons may change.
! On input, f(:,nlonp4,j0:j1,nflds) is defined for tasks with mytidi==0.
! On output, f(:,i0:i1,j0:j1,nflds) is defined for all tasks.
!
! Args:
      integer,intent(in) :: k0,k1,i0,i1,j0,j1,nflds
      real,intent(inout) :: f(k0:k1,nlonp4,j0:j1,nflds)
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i,k,j,n,nlevs,nlons,nlats,nlonsend,len,idest,isrc,ier,
     |  isend,irecv,itask,lonsend0,lonsend1,mtag
      real :: fmin,fmax
      real,allocatable,save ::
     |  sndbuf(:,:,:,:), ! send buffer (nlevs,mxlon,nlats,nflds)
     |  rcvbuf(:,:,:,:)  ! recv buffer (nlevs,mxlon,nlats,nflds)
#ifdef VT
!     call vtsymdef(107, 'mp_scatterlons_f3d','Communication',ier)
      call vtbegin(107,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_scatterlons_f3d,tsec,'begin')
!
! Allocate send and receive buffers:
      nlons = i1-i0+1
      nlevs = k1-k0+1
      nlats = j1-j0+1
      if (nflds==0) then
        write(6,"('>>> WARNING mp_scatterlons_f3d: no fields? nflds=',
     |    i3,' returning.')") nflds
        return
      endif
      len = nlevs*mxlon*mxlat*nflds
!
! 8/30/11 btf: If the buffers are allocated only once per run in 
!   this subroutine (i.e., below conditionals are opened up and 
!   deallocations at end are removed), then the results are different
!   (probably wrong) than when the buffers are allocated/deallocated 
!   in every call.
!
!     if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_scatterlons_f3d: error allocating sndbuf')")
        else
!         write(6,"('mp_scatterlons_f3d mytid=',i3,': allocated sndbuf',
!    |      ' nlevs=',i4,' mxlon=',i4,' mxlat=',i4,' nflds=',i4,' len=',
!    |      i8)") mytid,nlevs,mxlon,mxlat,nflds,len
        endif
!     endif
!     if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_scatterlons_f3d: error allocating rcvbuf')")
        else
!         write(6,"('mp_scatterlons_f3d mytid=',i3,': allocated rcvbuf',
!    |      ' nlevs=',i4,' mxlon=',i4,' mxlat=',i4,' nflds=',i4,' len=',
!    |      i8)") mytid,nlevs,mxlon,mxlat,nflds,len
        endif
!     endif
      sndbuf = 0.
      rcvbuf = 0.
!
! If mytidi==0, send to other tasks in my row (mytidi>0,mytidj):
      if (mytidi == 0) then
        do itask=1,ntaski-1
          idest = itask_table(itask,mytidj)
          lonsend0 = tasks(idest)%lon0
          lonsend1 = tasks(idest)%lon1
          nlonsend = lonsend1-lonsend0+1

!       write(6,"('mp_scatterlons mytidi=',i3,' idest=',i3,': i0,i1=',
!    |    2i4,' nlons=',i4,' j0,j1=',2i4,' nlats=',i4,' name=',a,
!    |    ' fin    min,max=',2e12.4)") 
!    |    mytidi,idest,lonsend0,lonsend1,nlonsend,j0,j1,nlats,name,
!    |    minval(f(:,lonsend0:lonsend1,j0:j1,1)),
!    |    maxval(f(:,lonsend0:lonsend1,j0:j1,1))

          mtag = idest+mytid
          do n=1,nflds
            do j=j0,j1
              sndbuf(:,1:nlonsend,j-j0+1,n) = f(:,lonsend0:lonsend1,j,n) 
            enddo ! j=j0,j1
          enddo ! n=1,nflds
          
!       write(6,"('mp_scatterlons mytidi=',i3,' idest=',i3,': i0,i1=',
!    |    2i4,' nlons=',i4,' j0,j1=',2i4,' nlats=',i4,' name=',a,
!    |    ' sndbuf min,max=',2e12.4)") 
!    |    mytidi,idest,lonsend0,lonsend1,nlonsend,j0,j1,nlats,name,
!    |    minval(sndbuf(:,1:nlonsend,1:nlats,1)),
!    |    maxval(sndbuf(:,1:nlonsend,1:nlats,1))

          mtag = idest+mytid
          call mpi_isend(sndbuf,len,MPI_REAL8,idest,mtag,MPI_COMM_WORLD,
     |      isend,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_scatterlons_f3d send to idest')
          call mpi_wait(isend,irstat,ier)
          if (ier /= 0) 
     |      call handle_mpi_err(ier,'mp_scatterlons_f3d wait for send')
        enddo ! itask=1,ntaski-1
!
! If mytidi > 0, receive from task (0,mytidj):
      else
        isrc = itask_table(0,mytidj)
        mtag = isrc+mytid
        call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,mtag,MPI_COMM_WORLD,
     |    irecv,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_scatterlons_f3d recv fm isrc')
        call mpi_wait(irecv,irstat,ier)
        if (ier /= 0) 
     |    call handle_mpi_err(ier,'mp_scatterlons_f3d wait for recv')
        do n=1,nflds
          do j=j0,j1 
            f(:,i0:i1,j,n) = rcvbuf(:,1:nlons,j-j0+1,n)
          enddo ! j=j0,j1
        enddo ! n=1,nflds
          
!       write(6,"('mp_scatterlons mytidi=',i3,' isrc =',i3,': i0,i1=',
!    |    2i4,' nlons=',i4,' j0,j1=',2i4,' nlats=',i4,' name=',a,
!    |    ' rcvbuf min,max=',2e12.4)") 
!    |    mytidi,isrc,i0,i1,nlons,j0,j1,nlats,name,
!    |    minval(rcvbuf(:,1:nlons,1:nlats,1)),
!    |    maxval(rcvbuf(:,1:nlons,1:nlats,1))
          
!       write(6,"('mp_scatterlons mytidi=',i3,' isrc =',i3,': i0,i1=',
!    |    2i4,' nlons=',i4,' j0,j1=',2i4,' nlats=',i4,' name=',a,
!    |    ' fout   min,max=',2e12.4)") 
!    |    mytidi,isrc,i0,i1,nlons,j0,j1,nlats,name,
!    |    minval(f(:,i0:i1,j0:j1,1)),maxval(f(:,i0:i1,j0:j1,1))

      endif
!
! Free local buffer space:
! 8/30/11 btf: See comments at allocation statements above.
      deallocate(sndbuf)
      deallocate(rcvbuf)
!
      if (do_rtc_mpi) then
        call timer(rtc0_scatterlons_f3d,tsec,'end')
        rtcmp_scatterlons_f3d = rtcmp_scatterlons_f3d+tsec
      endif
#ifdef VT
!     call vtsymdef(107, 'mp_scatterlons_f3d','Communication',ier)
      call vtend(107,ier)
#endif
      end subroutine mp_scatterlons_f3d
!-----------------------------------------------------------------------
      subroutine mp_periodic_f4d(ixt)
      use fields_module,only: i_tn ! for debug
!
! Define periodic points for all fields f4d(nf4d_hist)
!   (all latitudes and vertical column). This is called after
!   reading source history data.
! Periodic points:
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: idest,isrc,isend,irecv,ier,i,j,jj,n,len,k,nlats,mtag,
     |  status(MPI_STATUS_SIZE)
      real,allocatable,save :: sndbuf(:,:,:,:) ! (nlevp1,nlats,nf4d_hist,2)
      real,allocatable,save :: rcvbuf(:,:,:,:) ! (nlevp1,nlats,nf4d_hist,2)
!
#ifdef VT
!     call vtsymdef(108, 'mp_periodic_f4d','Communication',ier)
      call vtbegin(108,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_periodic_f4d,tsec,'begin')
!
! Use assignment and return in the special case of a serial run:
      if (ntask == 1) then     ! lon0,lon1 == 1,nlonp4
        do n=1,nf4d_hist       ! loop for fields on history at time ixt
          do j=lat0-2,lat1+2   ! latitude
            f4d(n)%data(:,lon0:lon0+1,j,ixt) =     ! e.g.: 1,2 <= 73,74
     |        f4d(n)%data(:,nlon+1:nlon+2,j,ixt)
            f4d(n)%data(:,lon1-1:lon1,j,ixt) =     ! e.g.: 75,76 <= 3,4
     |        f4d(n)%data(:,lon0+2:lon0+3,j,ixt)

!           write(6,"(/,'mp_periodic_f4d ntask=',i2,' j=',i3,' n=',i3,
!    |        ' field ',a,'(:,lon1-1,j)=',/,(6e12.4))") ntask,j,n,
!    |        trim(f4d(n)%short_name),f4d(n)%data(:,lon1-1,j,ixt)
!           write(6,"('  field ',a,'(:,lon1,j)=',/,(6e12.4))")
!    |        trim(f4d(n)%short_name),f4d(n)%data(:,lon1,j,ixt)

          enddo
        enddo
        return
      endif
!
! The restriction that every task must carry at least 4 longitudes
! means that tasks with mytidi==0 will have lons 1->4 and tasks 
! with mytidi==ntaski-1 will have lons nlonp4-3->nlonp4.
!
      if (mytidi /= 0 .and. mytidi /= ntaski-1) return
!
! Allocate send and receive buffers:
      nlats = mxlat+4 ! +4 is to cover lat0-1,lat0-2,lat1+1,lat1+2
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevp1,nlats,nf4d_hist,2),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f4d: error allocating sndbuf:',
     |      ' nlevp1=',i3,' nlats=',i3,' nf4d_hist=',i3,' ier=',i4)") 
     |      nlevp1,nlats,nf4d_hist,ier
      endif
!
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,nlats,nf4d_hist,2),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f4d: error allocating rcvbuf:',
     |      ' nlevp1=',i3,' nlats=',i3,' nf4d_hist=',i3,' ier=',i4)") 
     |      nlevp1,nlats,nf4d_hist,ier
      endif
      sndbuf = 0. ; rcvbuf = 0.
      len = nlevp1*nlats*nf4d_hist*2
!
! If I am at west end (mytidi==0) of task row, send lons 3,4 to task 
!   on east end (mytidi==ntaski-1) for its lons nlonp4-1,nlonp4.
!
      if (mytidi==0) then ! lon0==1, send 3,4
        idest = itask_table(ntaski-1,mytidj)
        do n=1,nf4d_hist
          jj = 1
          do j=lat0-2,lat1+2
            sndbuf(:,jj,n,1) = f4d(n)%data(:,lon0+2,j,ixt)
            sndbuf(:,jj,n,2) = f4d(n)%data(:,lon0+3,j,ixt)
            jj = jj+1
          enddo
        enddo
        isrc  = idest
        mtag = idest+mytid
        call mpi_sendrecv(sndbuf, len, MPI_REAL8, idest, mtag, rcvbuf,
     |       len, MPI_REAL8, isrc, isrc+mytid+1, MPI_COMM_WORLD, status,
     |       ier)
        if (ier /= 0) call handle_mpi_err(ier,
     |      'mp_periodic_f4d sendrecv from west end')
!
! If I am at east end of task row (mytidi==ntaski-1), send lons 
!   lon1-3,lon1-2 to task on west end (mytidi==0) for its lons 1,2.
!
      elseif (mytidi==ntaski-1) then ! lon1==nlonp4, send lon1-3,lon1-2
        idest = itask_table(0,mytidj)
        do n=1,nf4d_hist
          jj = 1
          do j=lat0-2,lat1+2
            sndbuf(:,jj,n,1) = f4d(n)%data(:,lon1-3,j,ixt)
            sndbuf(:,jj,n,2) = f4d(n)%data(:,lon1-2,j,ixt)
            jj = jj+1
          enddo
        enddo
        isrc  = idest
        mtag = idest+mytid+1
        call mpi_sendrecv(sndbuf, len, MPI_REAL8, idest, mtag, rcvbuf,
     |       len, MPI_REAL8, isrc, isrc+mytid, MPI_COMM_WORLD, status,
     |       ier)
        if (ier /= 0) call handle_mpi_err(ier,
     |      'mp_periodic_f4d sendrecv from east end')
      endif
!
! If I am on west end of task row, receive lons 1,2.
!
      if (mytidi==0) then ! lon0==1
        do n=1,nf4d_hist
          jj = 1
          do j=lat0-2,lat1+2
            f4d(n)%data(:,lon0  ,j,ixt) = rcvbuf(:,jj,n,1)
            f4d(n)%data(:,lon0+1,j,ixt) = rcvbuf(:,jj,n,2)
            jj = jj+1
          enddo
        enddo
!
! If I am on east end of task row, receive lons lon1-1,lon1:
!
      elseif (mytidi==ntaski-1) then ! lon1==nlonp4
        do n=1,nf4d_hist
          jj = 1
          do j=lat0-2,lat1+2
            f4d(n)%data(:,lon1-1,j,ixt) = rcvbuf(:,jj,n,1)
            f4d(n)%data(:,lon1  ,j,ixt) = rcvbuf(:,jj,n,2)
            jj = jj+1
          enddo
        enddo
      endif  
      if (do_rtc_mpi) then
        call timer(rtc0_periodic_f4d,tsec,'end')
        rtcmp_periodic_f4d = rtcmp_periodic_f4d+tsec
      endif
#ifdef VT
!     call vtsymdef(108, 'mp_periodic_f4d','Communication',ier)
      call vtend(108,ier)
#endif
      end subroutine mp_periodic_f4d
!-----------------------------------------------------------------------
      subroutine mp_periodic_f3d(f,lev0,lev1,lon0,lon1,lat0,lat1,nf)
!
! Define periodic points for field f(nlevels,lon0:lon1,lat0:lat1):
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,nf
      real,intent(inout) :: f(lev0:lev1,lon0:lon1,lat0:lat1,nf)
!
! Local:
      integer :: k,j,idest,isrc,isend,irecv,ier,len,nlevs,nlons,nlats,n
      real,allocatable,save :: sndbuf(:,:,:,:),rcvbuf(:,:,:,:) ! (nlevs,2,mxlat,nf)
      integer :: nfsave=0
!
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f3d','Communication',ier)
      call vtbegin(109,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_periodic_f3d,tsec,'begin')
#ifdef MPI
!
! Serial run of an MPI job (1 task):
      if (ntask==1) then
        do n=1,nf
          do j=lat0,lat1
            f(:,lon0:lon0+1,j,n) = f(:,nlon+1:nlon+2,j,n)
            f(:,lon1-1:lon1,j,n) = f(:,lon0+2:lon0+3,j,n)
          enddo ! j=lat0,lat1
        enddo
        return
      endif
#endif
!
! The restriction that every task must carry at least 4 longitudes
! means that tasks with mytidi==0 will have lons 1->4 and tasks 
! with mytidi==ntaski-1 will have lons nlonp4-3->nlonp4.
!
      if (mytidi /= 0 .and. mytidi /= ntaski-1) return
!
! Allocate send and receive buffers:
      nlevs = lev1-lev0+1
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
      if (nf > nfsave) then
!       write(6,"('mp_periodic_f3d: nf=',i3,' nfsave=',i3,' -- will ',
!    |    ' reallocate buffers..')") nf,nfsave
        if (allocated(sndbuf)) deallocate(sndbuf)
        if (allocated(rcvbuf)) deallocate(rcvbuf)
        nfsave = nf
      endif
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevs,2,mxlat,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f3d: error allocating sndbuf:',
     |      ' nlevs=',i3,' nlats=',i3,' nf=',i3,' ier=',i4)")
     |      nlevs,nlats,nf,ier
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevs,2,mxlat,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f3d: error allocating rcvbuf:',
     |      ' nlevs=',i3,' mxlat=',i3,' nf=',i3,' ier=',i4)")
     |      nlevs,mxlat,nf,ier
      endif
      sndbuf = 0. ; rcvbuf = 0.
      len = nlevs*2*mxlat*nf
!
! Send lons 3,4 to lons nlonp4-1,nlonp4 at task with mytidi==ntaski-1,
! and receive lons 1,2 from the same task.
      if (mytidi==0) then
        idest = itask_table(ntaski-1,mytidj)
        do n=1,nf
          do j=lat0,lat1
            do k=lev0,lev1
              sndbuf(k-lev0+1,1:2,j-lat0+1,n) = f(k,3:4,j,n)
            enddo
          enddo ! j=lat0,lat1
        enddo

        call mpi_isend(sndbuf(:,:,:,1:nf),len,MPI_REAL8,idest,1,
     |    MPI_COMM_WORLD,isend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f3d send to idest')
        isrc  = idest
        call mpi_irecv(rcvbuf(:,:,:,1:nf),len,MPI_REAL8,isrc,1,
     |    MPI_COMM_WORLD,irecv,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f3d recv fm isrc')
!
! Send lons nlonp4-3,nlonp4-2 to lons 1,2 at task with mytidi==0,
! and receive lons nlonp4-1,nlonp4 from same task:
      elseif (mytidi==ntaski-1) then
        idest = itask_table(0,mytidj)

        do n=1,nf
          do j=lat0,lat1
            do k=lev0,lev1
              sndbuf(k-lev0+1,1:2,j-lat0+1,n) =
     |          f(k,nlonp4-3:nlonp4-2,j,n)
            enddo
          enddo ! j=lat0,lat1
        enddo

        call mpi_isend(sndbuf(:,:,:,1:nf),len,MPI_REAL8,idest,1,
     |    MPI_COMM_WORLD,isend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f3d send to idest')
        isrc  = idest
        call mpi_irecv(rcvbuf(:,:,:,1:nf),len,MPI_REAL8,isrc,1,
     |    MPI_COMM_WORLD,irecv,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f3d recv fm isrc')
      endif
!
! Wait for completions:
      call mpi_wait(isend,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_periodic_f3d wait for send')
      call mpi_wait(irecv,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_periodic_f3d wait for recv')
!
! Copy from receive buffer:
      if (mytidi==0) then ! I am on west end
        do n=1,nf
          do j=lat0,lat1
            do k=lev0,lev1
              f(k,1:2,j,n) = rcvbuf(k-lev0+1,1:2,j-lat0+1,n)
            enddo ! k=lev0,lev1
          enddo ! j=lat0,lat1
        enddo
      elseif (mytidi==ntaski-1) then ! I am on east end
        do n=1,nf
          do j=lat0,lat1
            do k=lev0,lev1
              f(k,nlonp4-1:nlonp4,j,n) = rcvbuf(k-lev0+1,1:2,j-lat0+1,n)
            enddo ! k=lev0,lev1
          enddo ! j=lat0,lat1
        enddo
      endif
      if (do_rtc_mpi) then
        call timer(rtc0_periodic_f3d,tsec,'end')
        rtcmp_periodic_f3d = rtcmp_periodic_f3d+tsec
      endif
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f3d','Communication',ier)
      call vtend(109,ier)
#endif
      end subroutine mp_periodic_f3d
!-----------------------------------------------------------------------
      subroutine mp_periodic_f2d(f,lon0,lon1,lat0,lat1,nf)
!
! Define periodic points for field f(lon0:lon1,lat0:lat1):
! lons 1,2 <- nlonp4-3,nlonp4-2 and nlonp4-1,nlonp4 <- 3,4
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,nf
      real,intent(inout) :: f(lon0:lon1,lat0:lat1,nf)
!
! Local:
      integer :: j,idest,isrc,isend,irecv,ier,len,nlons,nlats,n
      real,allocatable,save :: sndbuf(:,:,:),rcvbuf(:,:,:) ! (2,nlats,nf)
      integer :: nfsave=0
!
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f2d','Communication',ier)
      call vtbegin(109,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_periodic_f2d,tsec,'begin')
!
! Use assignment and return in the special case of a serial run:
      if (ntask == 1) then
        do n=1,nf
          do j=lat0,lat1       ! latitude
            f(lon0:lon0+1,j,n) = f(nlon+1:nlon+2,j,n) ! e.g.: 1,2 <= 73,74
            f(lon1-1:lon1,j,n) = f(lon0+2:lon0+3,j,n) ! e.g.: 75,76 <= 3,4
          enddo
        enddo
        return
      endif
!
! The restriction that every task must carry at least 4 longitudes
! means that tasks with mytidi==0 will have lons 1->4 and tasks 
! with mytidi==ntaski-1 will have lons nlonp4-3->nlonp4.
!
      if (mytidi /= 0 .and. mytidi /= ntaski-1) return
!
! Allocate send and receive buffers:
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
      if (nf > nfsave) then
!       write(6,"('mp_periodic_f2d: nf=',i3,' nfsave=',i3,' -- will ',
!    |    ' reallocate buffers..')") nf,nfsave
        if (allocated(sndbuf)) deallocate(sndbuf)
        if (allocated(rcvbuf)) deallocate(rcvbuf)
        nfsave = nf
      endif
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(2,mxlat,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f2d: error allocating sndbuf:',
     |      ' nlats=',i3,' nf=',i3,' ier=',i4)") nlats,nf,ier
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(2,mxlat,nf),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_periodic_f2d: error allocating rcvbuf:',
     |      ' nlats=',i3,' nf=',i3,' ier=',i4)") nlats,nf,ier
      endif
      sndbuf = 0. ; rcvbuf = 0.
      len = 2*mxlat*nf
!
! Send lons 3,4 to lons nlonp4-1,nlonp4 at task with mytidi==ntaski-1,
! and receive lons 1,2 from the same task.
      if (mytidi==0) then
        idest = itask_table(ntaski-1,mytidj)
        do n=1,nf
          do j=lat0,lat1
            sndbuf(1:2,j-lat0+1,n) = f(3:4,j,n)
          enddo ! j=lat0,lat1
        enddo
        call mpi_isend(sndbuf(:,:,1:nf),len,MPI_REAL8,idest,1,
     |    MPI_COMM_WORLD,isend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f2d send to idest')
        isrc  = idest
        call mpi_irecv(rcvbuf(:,:,1:nf),len,MPI_REAL8,isrc,1,
     |    MPI_COMM_WORLD,irecv,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f2d recv fm isrc')
!
! Send lons nlonp4-3,nlonp4-2 to lons 1,2 at task with mytidi==0,
! and receive lons nlonp4-1,nlonp4 from same task:
      elseif (mytidi==ntaski-1) then
        idest = itask_table(0,mytidj)
        do n=1,nf
          do j=lat0,lat1
            sndbuf(1:2,j-lat0+1,n) = f(nlonp4-3:nlonp4-2,j,n)
          enddo ! j=lat0,lat1
        enddo
        call mpi_isend(sndbuf(:,:,1:nf),len,MPI_REAL8,idest,1,
     |    MPI_COMM_WORLD,isend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f2d send to idest')
        isrc  = idest
        call mpi_irecv(rcvbuf(:,:,1:nf),len,MPI_REAL8,isrc,1,
     |    MPI_COMM_WORLD,irecv,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f2d recv fm isrc')
      endif
!
! Wait for completions:
      call mpi_wait(isend,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_periodic_f2d wait for send')
      call mpi_wait(irecv,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_periodic_f2d wait for recv')
!
! Copy from receive buffer:
      if (mytidi==0) then
        do n=1,nf
          do j=lat0,lat1
            f(1:2,j,n) = rcvbuf(1:2,j-lat0+1,n)
          enddo ! j=lat0,lat1
        enddo
      elseif (mytidi==ntaski-1) then
        do n=1,nf
          do j=lat0,lat1
            f(nlonp4-1:nlonp4,j,n) = rcvbuf(1:2,j-lat0+1,n)
          enddo ! j=lat0,lat1
        enddo
      endif
      if (do_rtc_mpi) then
        call timer(rtc0_periodic_f2d,tsec,'end')
        rtcmp_periodic_f2d = rtcmp_periodic_f2d+tsec
      endif
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f2d','Communication',ier)
      call vtend(109,ier)
#endif
      end subroutine mp_periodic_f2d
!-----------------------------------------------------------------------
      subroutine mp_dynpot(ixt)
!
! Electric potential (geographic) is dimensioned in the fields module 
!   as dynpot(nlonp1,0:nlatp1,nlevp1). This routine is called from 
!   rdsource.F (after reading the source history) to define dynpot 
!   from the subdomains of f4d(poten) which were read from the source. 
!   Each task sends its f4d(poten) subdomain to all other tasks,
!   which fill in their dynpot.
!
! Args:
      integer,intent(in) :: ixt
!
! Local:
      integer :: i,j,ier,lonbeg,lonend,len,n,lat0n,lat1n
      real,allocatable,save :: sndbuf(:,:,:,:), rcvbuf(:,:,:,:)
      real :: fmin,fmax
#ifdef VT
!     call vtsymdef(110, 'mp_dynpot','Communication',ier)
      call vtbegin(110,ier)
#endif
      if (do_rtc_mpi) call timer(rtc0_dynpot,tsec,'begin')
!
! Allocate send and receive buffers for subdomains:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevp1,mxlon,mxlat,0:ntask-1),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_dynpot: error allocating sndbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' ier=',i4)") mxlon,mxlat,ier
          call shutdown('allocate in mp_dynpot')
        endif
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,mxlon,mxlat,0:ntask-1),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_dynpot: error allocating rcvbuf:',
     |      ' mxlon=',i3,' mxlat=',i3,' ier=',i4)") mxlon,mxlat,ier
          call shutdown('allocate in mp_dynpot')
        endif
      endif
      sndbuf = 0. ; rcvbuf = 0.
!
! Note len is the number of elements sent to *each* task, not the
!   total length of the snd/rcv buffers (the buffers contain data for 
!   *all* tasks). 
! In mpi_alltoall doc, "The jth block of data sent from task i is 
!   received by task j and placed in the ith block of the buffer recvbuf"
!
      len = nlevp1*mxlon*mxlat
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-1
!
! Load send buffer for all tasks with subdomain of current task:
      do n=0,ntask-1
        do j=lat0,lat1
          do i=lonbeg,lonend
            sndbuf(:,i-lonbeg+1,j-lat0+1,n) = 
     |        f4d(i_poten)%data(:,i,j,ixt)
          enddo ! i=lonbeg,lonend
        enddo ! j=lat0,lat1
      enddo ! n=0,ntask-1
!
! Do broadcast and receive in single alltoall mpi call:
!
      call mpi_alltoall(sndbuf,len,MPI_REAL8,
     |                  rcvbuf,len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_dynpot call mpi_alltoall')
!
! Define dynpot from receive buffers for each task (includes myself):
! lat0n,lat1n = latitudes of n'th subdomain,
! lonbeg,lonend = longitudes of n'th subdomain 
!   (excluding end points 1,2 and nlonp4-1,nlonp4-2)
!
      do n=0,ntask-1
        lat0n = tasks(n)%lat0
        lat1n = tasks(n)%lat1
        lonbeg = tasks(n)%lon0
        if (lonbeg==1) lonbeg = 3
        lonend = tasks(n)%lon1
        if (lonend==nlonp4) lonend = nlonp4-1
!
! Loop over n'th subdomain:
        do j=lat0n,lat1n
          do i=lonbeg,lonend
            dynpot(i-2,j,:) = rcvbuf(:,i-lonbeg+1,j-lat0n+1,n)
          enddo ! i=lonbeg,lonend
        enddo ! j=lat0,lat1
      enddo ! n=0,ntask-1
      if (do_rtc_mpi) then
        call timer(rtc0_dynpot,tsec,'end')
        rtcmp_dynpot = rtcmp_dynpot+tsec
      endif
#ifdef VT
!     call vtsymdef(110, 'mp_dynpot','Communication',ier)
      call vtend(110,ier)
#endif
      end subroutine mp_dynpot
!-----------------------------------------------------------------------
      subroutine mp_updatephi
!
! Dynamo has been calculated by the master task. Now send dynamo output
!   field phim3d (mag coords) to the non-root tasks . Phim3d is
!   dimensioned at the global mag grid by all tasks (see fields.F):
!   phim3d(nmlonp1,nmlat,nmlev), but it is calculated in the dynamo by
!   only the master task.
!
! 8/31/11 btf (as in timegcm 11/4/09): Using a simple mpi_bcast call here 
!   instead of mpi_send/mpi_recv is *much* faster, e.g., with 64 procs 
!   running timegcm on bluefire, bcast showed 12% speedup overall at 
!   double-res, and 25% speedup at single-res.
!
! Local:
      integer :: ier,len
      real :: fmin,fmax
!
      len = nmlonp1*nmlat*nmlev
      call mpi_bcast(phim3d,len,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_updatephi call mpi_bcast')
!     if (mytid /= 0) then
!       call fminmax(phim3d,len,fmin,fmax)
!       write(6,"('mp_updatephi: task ',i3,' received phim3d ',
!    |    'from master: global min,max=',2e12.4)") mytid,fmin,fmax
!     endif
      end subroutine mp_updatephi      
!-----------------------------------------------------------------------
      subroutine mp_updateemphi
!
! Dynamo has been calculated by the master task. Now send dynamo output
!   field emphi3d (mag coords) to the slaves (emphi3d is use-associated
!   from fields module at top of this module). Emphi3d is dimensioned at
!   the global mag grid by all tasks: emphi3d(nmlonp1,nmlat,nlevp1),
!   but it is calculated in the dynamo by only the master task.  
!
! Local:
      integer :: n,msgtag,ier
      real :: fmin,fmax
!
      if (mytid==0) then ! master send to slaves
        do n=1,ntask-1
          msgtag = n
          call mpi_send(emphi3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,n,
     |      msgtag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
            write(6,"('>>> mp_updatephi: error sending',
     |        ' emphi3d from master to task ',i3,' msgtag=',i3)")
     |        n,msgtag
          else
!           write(6,"('mp_updatephi: master sent emphi3d to task',
!    |        i3)") n
          endif
        enddo
      else               ! slave receive from master
        msgtag = mytid
        call mpi_recv(emphi3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,0,
     |    msgtag,MPI_COMM_WORLD,irstat,ier)
        if (ier /= 0) then
          write(6,"('>>> mp_updatephi: error receiving',
     |      ' emphi3d from master at task mytid=',i3,' msgtag=',i3)")
     |      mytid,msgtag
        endif
      endif
      end subroutine mp_updateemphi
!-----------------------------------------------------------------------
      subroutine mp_updateemlam
!
! Dynamo has been calculated by the master task. Now send dynamo output
!   field emlam3d (mag coords) to the slaves (emlam3d is use-associated
!   from fields module at top of this module). Emlam3d is dimensioned at
!   the global mag grid by all tasks: emlam3d(nmlonp1,nmlat,nlevp1),
!   but it is calculated in the dynamo by only the master task.  
!
! Local:
      integer :: n,msgtag,ier
      real :: fmin,fmax
!
      if (mytid==0) then ! master send to slaves
        do n=1,ntask-1
          msgtag = n
          call mpi_send(emlam3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,n,
     |      msgtag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
            write(6,"('>>> mp_updatemlam: error sending',
     |        ' emlam3d from master to task ',i3,' msgtag=',i3)")
     |        n,msgtag
          else
!           write(6,"('mp_updatemlam: master sent emlam3d to task',
!    |        i3)") n
          endif
        enddo
      else               ! slave receive from master
        msgtag = mytid
        call mpi_recv(emlam3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,0,
     |    msgtag,MPI_COMM_WORLD,irstat,ier)
        if (ier /= 0) then
          write(6,"('>>> mp_updatemlam: error receiving',
     |      ' emlam3d from master at task mytid=',i3,' msgtag=',i3)")
     |      mytid,msgtag
        endif
      endif
      end subroutine mp_updateemlam      
!-----------------------------------------------------------------------
      subroutine mp_updateemz
!
! Dynamo has been calculated by the master task. Now send dynamo output
!   field emz3d (mag coords) to the slaves (emz3d is use-associated
!   from fields module at top of this module). Emz3d is dimensioned at
!   the global mag grid by all tasks: emz3d(nmlonp1,nmlat,nlevp1),
!   but it is calculated in the dynamo by only the master task.  
!
! Local:
      integer :: n,msgtag,ier
      real :: fmin,fmax
!
      if (mytid==0) then ! master send to slaves
        do n=1,ntask-1
          msgtag = n
          call mpi_send(emz3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,n,
     |      msgtag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
            write(6,"('>>> mp_updateemz: error sending',
     |        ' emz3d from master to task ',i3,' msgtag=',i3)")
     |        n,msgtag
          else
!           write(6,"('mp_updatemz: master sent emz3d to task',
!    |        i3)") n
          endif
        enddo
      else               ! slave receive from master
        msgtag = mytid
        call mpi_recv(emz3d,nmlonp1*nmlat*(nlevp1+3),MPI_REAL8,0,
     |    msgtag,MPI_COMM_WORLD,irstat,ier)
        if (ier /= 0) then
          write(6,"('>>> mp_updateemz: error receiving',
     |      ' emz3d from master at task mytid=',i3,' msgtag=',i3)")
     |      mytid,msgtag
        endif
      endif
      end subroutine mp_updateemz
!-----------------------------------------------------------------------
      subroutine mp_close
      integer :: ier
!
      call mpi_finalize(ier)
      if (ier /= 0) then
        write(6,"(/,'>>> WARNING: error from mp_finalize: ier=',i3)")ier
      endif
      end subroutine mp_close
!-----------------------------------------------------------------------
      subroutine handle_mpi_err(ierrcode,string)
!
! Args:
      integer,intent(in) :: ierrcode
      character(len=*) :: string
!
! Local:
      character(len=80) :: errstring
      integer :: len_errstring
!
      call mpi_error_string(ierrcode,errstring,len_errstring)
      write(6,"(/,'>>> mpi error: ',a)") trim(string)
      write(6,"('  ierrcode=',i3,': ',a)") trim(errstring)
      end subroutine handle_mpi_err
#endif
      end module mpi_module
