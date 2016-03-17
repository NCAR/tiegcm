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
     |  nf4d,nf4d_hist,fsechist,foutput,tlbc,ulbc,vlbc,tlbc_glb,
     |  ulbc_glb,vlbc_glb,emphi3d,emlam3d,emz3d,fzg,tlbc_nm,ulbc_nm,
     |  vlbc_nm,tlbc_nm_glb,ulbc_nm_glb,vlbc_nm_glb
#ifdef MPI
      use mpi
#endif
      use hist_module,only: nfsech
      use input_module,only: ntask_lat,ntask_lon,
     |  ntask_maglat,ntask_maglon
      implicit none
#ifdef MPI
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
! Geographic subdomains for current task:
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
     |  itask_table(:,:)  ! 2d table of tasks (i,j)
!
! Magnetic subdomains for current task:
      integer ::
     |  nmagtaski,        ! number of tasks in mag lon dimension
     |  nmagtaskj,        ! number of tasks in mag lat dimension
     |  magtidi,          ! i coord for current task in task table
     |  magtidj,          ! j coord for current task in task table
     |  mlat0,mlat1,      ! first and last lats for each task
     |  mlon0,mlon1,      ! first and last lons for each task
     |  mxmaglon,mxmaglat ! max number of mag subdomain lon,lat points among all tasks
      integer,allocatable,save ::
     |  itask_table_mag(:,:)  ! 2d table of tasks (i,j)
      integer :: cols_comm  ! communicators for each task column
!
! Task type:
      type task
        integer :: mytid     ! task id
!
! Geographic subdomains:
        integer :: mytidi    ! task coord in longitude dimension of task table
        integer :: mytidj    ! task coord in latitude  dimension of task table
        integer :: nlats     ! number of latitudes  calculated by this task
        integer :: nlons     ! number of longitudes calculated by this task
        integer :: lat0,lat1 ! first and last latitude  indices
        integer :: lon0,lon1 ! first and last longitude indices
!
! Magnetic subdomains:
        integer :: magtidi        ! task coord in mag longitude dimension of task table
        integer :: magtidj        ! task coord in mag latitude  dimension of task table
        integer :: nmaglats       ! number of mag latitudes  calculated by this task
        integer :: nmaglons       ! number of mag longitudes calculated by this task
        integer :: mlat0,mlat1    ! first and last latitude  indices
        integer :: mlon0,mlon1    ! first and last longitude indices
      end type task
!
! Conjugate points in mag subdomains, for mp_mag_foldhem
      integer,allocatable,dimension(:),save ::    ! (ntask)
     |  nsend_south,       ! number of south lats to send to north (each task)
     |  nrecv_north        ! number of north lats to send to south (each task)
      integer,allocatable,dimension(:,:),save ::  ! (mxlats,ntask)
     |  send_south_coords, ! south j lats to send to north
     |  recv_north_coords  ! north j lats to recv from south
!
! type(task) :: tasks(ntask) will be made available to all tasks
! (so each task has information about all tasks)
!
      type(task),allocatable :: tasks(:)
!
! MPI timing:
      logical :: mpi_timing=.true.
      real,save ::
     |  time_totalrun = 0.,
     |  time_gather2root_prim = 0.,
     |  time_gather2root_f3d = 0.,
     |  time_gather2root_sech = 0.,
     |  time_gather2root_lbc = 0.,
     |  time_bndlats = 0.,
     |  time_bndlats_f2d = 0.,
     |  time_bndlons = 0.,
     |  time_bndlons_f3d = 0.,
     |  time_bndlons_f2d = 0.,
     |  time_polelats = 0.,
     |  time_polelat_f3d = 0.,
     |  time_gatherlons_f3d = 0.,
     |  time_scatterlons_f3d = 0.,
     |  time_periodic_f4d = 0.,
     |  time_periodic_f3d = 0.,
     |  time_periodic_f2d = 0.,
     |  time_mageq = 0.,
     |  time_mageq_jpm1 = 0.,
     |  time_mageq_jpm3 = 0.,
     |  time_magpole_2d = 0.,
     |  time_magpole_3d = 0.,
     |  time_magpoles = 0.,
     |  time_conjugate_points = 0.,
     |  time_mag_foldhem = 0.,
     |  time_mag_periodic_f2d = 0.,
     |  time_gather_pdyn = 0.,
     |  time_mag_halos = 0.,
     |  time_geo_halos = 0.,
     |  time_geo_halos_f3d = 0.,
     |  time_scatter_coeffs = 0.,
     |  time_scatter_phim = 0.,
     |  time_bndlats_kmh = 0.,
     |  time_bndlons_kmh = 0.,
     |  time_mag_jslot = 0.,
     |  time_gather_f2d = 0.,
     |  time_scatter_f2d = 0.
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
      integer,dimension(11) :: 
     |  valid_pecounts=(/1,4,8,12,16,24,32,48,64,72,80/)
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
     |    'has been validated:',/,(12i6))") ntask,valid_pecounts
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
!     write(6,"('mp_init: ntask=',i3,' mytid=',i3)") ntask,mytid
      end subroutine mp_init
!-----------------------------------------------------------------------
      subroutine mp_distribute_geo
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
!
!   tasks(ntasks): array of task types.
!   (note mytidi,j,lat0,1,lon0,1, above are redundant with
!    tasks(mytid), eg, tasks(mytid)%mytidi==mytidi)
!
! Local:
      integer :: i,j,n,ier,irank,nj,ni,ncells
!
! Set ntaski,j from input ntask_lon,ntask_lat:
      ntaski = ntask_lon
      ntaskj = ntask_lat
!
! ntaski*ntaskj must == ntask
! (This was already checked in input, but you never know...)
      if (ntaski*ntaskj /= ntask) then 
        write(6,"('>>> mp_distribute_geo: ntaski*ntaskj must == ',
     |    'ntask.')")
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
!
! Tasks are periodic in longitude:
! 8/5/14 btf: this was apparently causing a problem at lon=0,
!   esp visible in maps of UN and OMEGA at higher levels, but 
!   only when ntask==4.
! 
!       itask_table(-1,j) = itask_table(ntaski-1,j)
!       itask_table(ntaski,j) = itask_table(0,j)
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
      enddo
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
      end subroutine mp_distribute_geo
!-----------------------------------------------------------------------
      subroutine mp_distribute_mag
!
! Set up 2-d data decomposition for geomagnetic grid
! (see also sub mp_distribute_geo)
!
! Local:
      integer :: i,j,n,ier,irank,nj,ni,ncells,tidcol
!
! Set nmagtaski,j from input ntask_maglon,ntask_maglat:
      nmagtaski = ntask_maglon
      nmagtaskj = ntask_maglat
!
! nmagtaski*nmagtaskj must == ntask
! (This was already checked in input, but you never know...)
      if (nmagtaski*nmagtaskj /= ntask) then
        write(6,"('>>> mp_distribute_mag: nmagtaski*nmagtaskj must == ',
     |    'ntask.')")
        write(6,"('    nmagtaski=',i3,' nmagtaskj=',i3,' ntask=',i3)")
     |    nmagtaski,nmagtaskj,ntask
        call shutdown('nmagtaski*j')
      endif
!
! Allocate and set 2d table of tasks:
      allocate(itask_table_mag(-1:nmagtaski+1,-1:nmagtaskj+1),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> Error allocating itable: nmagtaski,j=',2i3)")
     |    nmagtaski,nmagtaskj
        call mp_close
        call shutdown('itask_table')
      endif
      itask_table_mag(:,:) = MPI_PROC_NULL
      irank = 0
      do j =  0,nmagtaskj-1
        do i = 0,nmagtaski-1
          itask_table_mag(i,j) = irank
          if (mytid == irank) then
            magtidi = i
            magtidj = j
          endif
          irank = irank+1
        enddo
!
! Tasks are periodic in longitude:
        itask_table_mag(-1,j) = itask_table_mag(nmagtaski-1,j)
        itask_table_mag(nmagtaski,j) = itask_table_mag(0,j)
      enddo
!
! Print table to stdout:
      write(6,"(/,'ntask=',i3,' ntask_maglon=',i2,' ntask_maglat=',i2,
     |  ' Task Table:')") ntask,ntask_maglon,ntask_maglat
      do j=-1,nmagtaskj
        write(6,"('j=',i3,' itask_table_mag(:,j)=',100i3)")
     |      j,itask_table_mag(:,j)
      enddo
!
! Calculate start and end indices in lon,lat dimensions for each task:
      call distribute_1d(1,nmlonp1,nmagtaski,magtidi,mlon0,mlon1)
      call distribute_1d(1,nmlat  ,nmagtaskj,magtidj,mlat0,mlat1)
      nj = mlat1-mlat0+1 ! number of mag latitudes for this task
      ni = mlon1-mlon0+1 ! number of mag longitudes for this task
      ncells = nj*ni   ! total number of grid cells for this task
!
! Report my stats to stdout:
      write(6,"(/,'mytid=',i3,' magtidi,j=',2i3,' mlat0,1=',2i3,' (',
     |  i2,') mlon0,1=',2i3,' (',i2,') ncells=',i4)") mytid,magtidi,
     |  magtidj,mlat0,mlat1,nj,mlon0,mlon1,ni,ncells
!
! Define all task structures with current task values
! (redundant for alltoall):
!
      do n=0,ntask-1
!       tasks(n)%mytid  = mytid      ! was set in mp_distribute_geo
        tasks(n)%magtidi = magtidi
        tasks(n)%magtidj = magtidj
        tasks(n)%nmaglats  = nj
        tasks(n)%nmaglons  = ni
        tasks(n)%mlat0   = mlat0
        tasks(n)%mlat1   = mlat1
        tasks(n)%mlon0   = mlon0
        tasks(n)%mlon1   = mlon1
      enddo
!
! All tasks must have at least 4 longitudes:
      do n=0,ntask-1
        if (tasks(n)%nmaglons < 4) then
          write(6,"(/,'>>> mp_distribute_mag: each task must carry ',
     |      'at least 4 longitudes. task=',i2,' nmaglons=',i2)")
     |      n,tasks(n)%nmaglons
          call shutdown('nmaglons per task')
        endif
      enddo
!
! Create subgroup communicators for each task column:
! These communicators will be used by sub mp_mag_jslot (mpi.F).
      call mpi_comm_split(MPI_COMM_WORLD,mod(mytid,nmagtaski),mytid,
     |  cols_comm,ier)
      call MPI_Comm_rank(cols_comm,tidcol,ier)

!     write(6,"('mp_distribute_geo: cols_comm=',i3,' nmagtaski=',i3,
!    |  ' tidcol=',i3)") cols_comm,nmagtaski,tidcol
      write(6,"('mp_distribute_geo: nmagtaski=',i3,' tidcol=',i3)")
     |  nmagtaski,tidcol

      end subroutine mp_distribute_mag
!-----------------------------------------------------------------------
      subroutine mp_exchange_tasks(iprint)
!
! Args:
      integer,intent(in) :: iprint
!
! Local:
! itasks_send(len_task_type,ntask) will be used to send tasks(:) info
!   to all tasks (directly passing mpi derived data types is reportedly
!   not stable, or not available until MPI 2.x).
!
      integer :: n,ier,i
      integer,parameter :: len_task_type = 17 ! see type task above
      integer,allocatable,save ::
     |  itasks_send(:,:), ! send buffer
     |  itasks_recv(:,:)  ! send buffer
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

        itasks_send(10,n) = tasks(mytid)%magtidi
        itasks_send(11,n) = tasks(mytid)%magtidj
        itasks_send(12,n) = tasks(mytid)%nmaglats
        itasks_send(13,n) = tasks(mytid)%nmaglons
        itasks_send(14,n) = tasks(mytid)%mlat0
        itasks_send(15,n) = tasks(mytid)%mlat1
        itasks_send(16,n) = tasks(mytid)%mlon0
        itasks_send(17,n) = tasks(mytid)%mlon1
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

        tasks(n)%magtidi   = itasks_recv(10,n)
        tasks(n)%magtidj   = itasks_recv(11,n)
        tasks(n)%nmaglats  = itasks_recv(12,n)
        tasks(n)%nmaglons  = itasks_recv(13,n)
        tasks(n)%mlat0   = itasks_recv(14,n)
        tasks(n)%mlat1   = itasks_recv(15,n)
        tasks(n)%mlon0   = itasks_recv(16,n)
        tasks(n)%mlon1   = itasks_recv(17,n)

!
        if (n==mytid.and.iprint > 0) then
          write(6,"(/,'Task ',i3,':')") n
          write(6,"(/,'Subdomain on geographic grid:')")
          write(6,"('tasks(',i3,')%mytid =',i3)") n,tasks(n)%mytid
          write(6,"('tasks(',i3,')%mytidi=',i3)") n,tasks(n)%mytidi
          write(6,"('tasks(',i3,')%mytidj=',i3)") n,tasks(n)%mytidj
          write(6,"('tasks(',i3,')%nlats =',i3)") n,tasks(n)%nlats
          write(6,"('tasks(',i3,')%nlons =',i3)") n,tasks(n)%nlons
          write(6,"('tasks(',i3,')%lat0  =',i3)") n,tasks(n)%lat0
          write(6,"('tasks(',i3,')%lat1  =',i3)") n,tasks(n)%lat1
          write(6,"('tasks(',i3,')%lon0  =',i3)") n,tasks(n)%lon0
          write(6,"('tasks(',i3,')%lon1  =',i3)") n,tasks(n)%lon1
          write(6,"('Number of geo subdomain grid points = ',i6)")
     |      tasks(n)%nlons * tasks(n)%nlats
          write(6,"(/,'Subdomain on geomagnetic grid:')")
          write(6,"('tasks(',i3,')%magtidi=',i3)") n,tasks(n)%magtidi
          write(6,"('tasks(',i3,')%magtidj=',i3)") n,tasks(n)%magtidj
          write(6,"('tasks(',i3,')%nmaglats =',i3)") n,tasks(n)%nmaglats
          write(6,"('tasks(',i3,')%nmaglons =',i3)") n,tasks(n)%nmaglons
          write(6,"('tasks(',i3,')%mlat0  =',i3)") n,tasks(n)%mlat0
          write(6,"('tasks(',i3,')%mlat1  =',i3)") n,tasks(n)%mlat1
          write(6,"('tasks(',i3,')%mlon0  =',i3)") n,tasks(n)%mlon0
          write(6,"('tasks(',i3,')%mlon1  =',i3)") n,tasks(n)%mlon1
          write(6,"('Number of mag subdomain grid points = ',i6)")
     |      tasks(n)%nmaglons * tasks(n)%nmaglats

!         write(6,"('Range of mag lats (deg) = ',f10.3,' Mag lats=',/,
!    |      (6f12.3))") gmlat(tasks(n)%mlat1)-gmlat(tasks(n)%mlat0),
!    |      gmlat(tasks(n)%mlat0:tasks(n)%mlat1)
!         write(6,"('Range of mag lons (deg) = ',f10.3,' Mag lons=',/,
!    |      (6f12.3))") gmlon(tasks(n)%mlon1)-gmlon(tasks(n)%mlon0),
!    |      gmlon(tasks(n)%mlon0:tasks(n)%mlon1)
        endif
      enddo
!
! Release locally allocated space:
! debug: temporarilly remove these deallocs while using getoverstepped()
      deallocate(itasks_send)
      deallocate(itasks_recv)
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
!     write(6,"('mp_exchange_tasks: mxlat=',i4,' mxlon=',i4)")
!    |  mxlat,mxlon
!
! mxmaglon,mxmaglat are maximum number of mag lons,lats owned by all tasks:
      mxmaglon = -9999
      do n=0,ntask-1
        if (tasks(n)%nmaglons > mxmaglon) mxmaglon = tasks(n)%nmaglons
      enddo
      mxmaglat = -9999
      do n=0,ntask-1
        if (tasks(n)%nmaglats > mxmaglat) mxmaglat = tasks(n)%nmaglats
      enddo
      write(6,"('mp_exchange_tasks: mxmaglat=',i4,' mxmaglon=',i4)")
     |  mxmaglat,mxmaglon
!
! Find conjugate points for folding hemispheres:
      call conjugate_points

      end subroutine mp_exchange_tasks
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
! Gather for either primary and secondary histories.
! If primary, also gather f3d field subdomains zg to global fzg.
      if (trim(type)=='prim') then
         call mp_gather2root_prim(ixt)
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
      real :: starttime,endtime
!
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gather2root_prim=time_gather2root_prim+(endtime-starttime)
      endif
      end subroutine mp_gather2root_prim
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
     |  ilat0,ilat1,msgtag,nlons,nlats,nlevs,maxlon,maxlat,maxlev
      integer :: idest = 0
      real :: fmin,fmax
      real :: starttime,endtime
      real,allocatable,save :: rcvbuf(:,:,:,:),sndbuf(:,:,:,:)

      if (mpi_timing) starttime = mpi_wtime()
!      
! Allocate send and receive buffers:
! Secondary history data fsechhist(i)%data(lon,lat,lev) (mag or geo)
!
      maxlon = max(mxlon,mxmaglon)
      maxlat = max(mxlat,mxmaglat)
      maxlev = max(nlevp1,nmlevp1)
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(maxlon,maxlat,maxlev,nfsech),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_sech: error allocating rcvbuf:',
     |      ' maxlon=',i3,' maxlat=',i3,' maxlev=',i3,' nfsech=',
     |      i3,' ier=',i4)") maxlon,maxlat,maxlev,nfsech,ier
      endif
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(maxlon,maxlat,maxlev,nfsech),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_gather2root_sech: error allocating sndbuf:',
     |      ' maxlon=',i3,' maxlat=',i3,' maxlev=',i3,' nfsech=',
     |      i3,' ier=',i4)") maxlon,maxlat,maxlev,nfsech,ier
      endif
      rcvbuf = 0. ; sndbuf = 0.
      len = maxlon*maxlat*maxlev*nfsech ! buffer length
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
! Load secondary history fields from receive buffer:
          do i = 1,nfsech
            if (.not.fsechist(i)%prognostic.and. 
     |          .not.fsechist(i)%task0_only.and.
     |          associated(fsechist(i)%data)) then
!
! Get bounds according to geo or mag:
              if (fsechist(i)%geo) then
                ilon0 = tasks(n)%lon0
                ilon1 = tasks(n)%lon1
                ilat0 = tasks(n)%lat0
                ilat1 = tasks(n)%lat1
                nlevs = nlevp1
              else
                ilon0 = tasks(n)%mlon0
                ilon1 = tasks(n)%mlon1
                ilat0 = tasks(n)%mlat0
                ilat1 = tasks(n)%mlat1
                nlevs = nmlevp1
              endif
              nlons = ilon1-ilon0+1
              nlats = ilat1-ilat0+1
!
! 3d field (lon,lat,lev):
              if (fsechist(i)%ndims==3) then
                fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,1:nlevs)=
     |                    rcvbuf(1:nlons,1:nlats,1:nlevs,i)
              else
!
! 2d field (lon,lat):
                fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,1)=
     |                    rcvbuf(1:nlons,1:nlats,1,i)
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
        do i = 1,nfsech
          if (.not.fsechist(i)%prognostic.and.
     |        .not.fsechist(i)%task0_only.and.
     |        associated(fsechist(i)%data)) then
            if (fsechist(i)%geo) then
              ilon0 = lon0
              ilon1 = lon1
              ilat0 = lat0
              ilat1 = lat1
              nlevs = nlevp1
            else
              ilon0 = mlon0
              ilon1 = mlon1
              ilat0 = mlat0
              ilat1 = mlat1
              nlevs = nmlevp1
            endif
            nlons = ilon1-ilon0+1
            nlats = ilat1-ilat0+1
!
! 3d send buffer:
            if (fsechist(i)%ndims==3) then     ! 3d send buffer
              sndbuf(1:nlons,1:nlats,1:nlevs,i) =
     |          fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,1:nlevs)
            else                               ! 2d send buffer
!
! 2d send buffer:
              sndbuf(1:nlons,1:nlats,1,i) =
     |          fsechist(i)%data(ilon0:ilon1,ilat0:ilat1,1)
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gather2root_sech=time_gather2root_sech+(endtime-starttime)
      endif
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
      real :: starttime,endtime
      real,allocatable,save :: rcvbuf(:,:,:),sndbuf(:,:,:)

      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gather2root_lbc=time_gather2root_lbc+(endtime-starttime)
      endif
      end subroutine mp_gather2root_lbc
!-----------------------------------------------------------------------
      subroutine mp_gather_f2d(f2d_sub,f2d_glb,lb1,ub1,lb2,ub2,
     |  ntasks,nglb1,nglb2,nf)
!
! Gather 2d array to root task. 
! On input, f2d_sub is 2d array on task subdomains
! On output, f3d_glb is 2d global array, defined only on the root task
!
! Called by sub laplacian prior to fft for HE upper boundary:
!     call mp_gather_f2d(lap_arg,fx,lon0,lon1,lat0,lat1,
!    |  nlonp4,nlat)
!
! Args:
      integer,intent(in) ::
     |  ntasks,        ! size of comm_world (should be = ntask)
     |  nglb1,nglb2,   ! size of global 1st and 2nd dimensions
     |  nf             ! number of arrays/fields
      integer,dimension(0:ntasks-1),intent(in) ::
     |  lb1,ub1,       ! lower,upper bound of 1st dimension of subdomains
     |  lb2,ub2        ! lower,upper bound of 2nd dimension of subdomains
! Input data at current subdomain:
      real,intent(in)  :: 
     |  f2d_sub(lb1(mytid):ub1(mytid),lb2(mytid):ub2(mytid),nf) 
      real,intent(out) :: f2d_glb(nglb1,nglb2,nf)     ! global output at root task
!
! Local:
      integer :: len,root,n,i,j,id,mxsub1,mxsub2,ier
      real,allocatable,save :: sndbuf(:,:,:),rcvbuf(:,:,:,:)
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()

      root = 0               ! gather to root task
      f2d_glb = 0.           ! init output
!
! Find max dimensions of all subdomains (might be mxlon,mxlat):
! (this only needs to be done once, but do it always for now)
!
      mxsub1 = 0. ; mxsub2 = 0.
      do id=0,ntasks-1
        len = ub1(id)-lb1(id)+1
        if (len > mxsub1) mxsub1 = len
        len = ub2(id)-lb2(id)+1
        if (len > mxsub2) mxsub2 = len
      enddo

      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxsub1,mxsub2,nf))
        sndbuf = 0.
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxsub1,mxsub2,nf,0:ntasks-1))
        rcvbuf = 0.
      endif
!
! Load send buffer with current subdomain:
      do n=1,nf
        do j=lb2(mytid),ub2(mytid)
          do i=lb1(mytid),ub1(mytid)
            sndbuf(i-lb1(mytid)+1,j-lb2(mytid)+1,n) = f2d_sub(i,j,n)
          enddo
        enddo
      enddo
!
! Do the gather:      
      len = mxsub1*mxsub2*nf ! max size of 2d subdomain
      call mpi_gather(sndbuf,len,MPI_REAL8,rcvbuf,len,MPI_REAL8,root,
     |  MPI_COMM_WORLD,ier)
!
! Define global output on the root task from receive buffer of each task:
      if (mytid==root) then
        do id=0,ntasks-1
          do n=1,nf
            do j=lb2(id),ub2(id)
              do i=lb1(id),ub1(id)
                f2d_glb(i,j,n) = rcvbuf(i-lb1(id)+1,j-lb2(id)+1,n,id)
              enddo
            enddo 
          enddo
        enddo
      endif      
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gather_f2d=time_gather_f2d+(endtime-starttime)
      endif
      end subroutine mp_gather_f2d
!-----------------------------------------------------------------------
      subroutine mp_scatter_f2d(f2d_glb,f2d_sub,lb1,ub1,lb2,ub2,
     |  ntasks,nglb1,nglb2,nf)
!
! Scatter global arrays on the root task to all task subdomains.
! This is essentially the reverse of mp_gather_f2d.
!
! Called from laplacian after fft for HE upper boundary:
!     call mp_scatter_f2d(fx,lap_out,tasks(:)%lon0,tasks(:)%lon1,
!    |  tasks(:)%lat0,tasks(:)%lat1,ntask,nlonp4,nlat,1)
!
! Args:
      integer,intent(in) ::
     |  ntasks,        ! size of comm_world (should be = ntask)
     |  nglb1,nglb2,   ! size of global 1st and 2nd dimensions
     |  nf             ! number of arrays/fields
      integer,dimension(0:ntasks-1),intent(in) ::
     |  lb1,ub1,       ! lower,upper bound of 1st dimension of subdomains
     |  lb2,ub2        ! lower,upper bound of 2nd dimension of subdomains
!
! f2d_glb is global input at root task
! f2d_sub is subdomain output on all tasks
!
      real,intent(in) :: f2d_glb(nglb1,nglb2,nf) 
      real,intent(out)  :: 
     |  f2d_sub(lb1(mytid):ub1(mytid),lb2(mytid):ub2(mytid),nf) 
!
! Local:
      integer :: len,root,n,i,j,id,mxsub1,mxsub2,ier
      real,allocatable,save :: sndbuf(:,:,:,:),rcvbuf(:,:,:)
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()

      root = 0     ! scatter from root task
      f2d_sub = 0. ! init output
!
! Find max dimensions of all subdomains (might be mxlon,mxlat):
! (this only needs to be done once, but do it always for now)
!
      mxsub1 = 0. ; mxsub2 = 0.
      do id=0,ntasks-1
        len = ub1(id)-lb1(id)+1
        if (len > mxsub1) mxsub1 = len
        len = ub2(id)-lb2(id)+1
        if (len > mxsub2) mxsub2 = len
      enddo

      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxsub1,mxsub2,nf,0:ntasks-1))
        sndbuf = 0.
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxsub1,mxsub2,nf))
        rcvbuf = 0.
      endif
!
! Load global send buffer from root task:
      if (mytid==root) then
        do id=0,ntasks-1
          do n=1,nf
            do j=lb2(id),ub2(id)
              do i=lb1(id),ub1(id)
                sndbuf(i-lb1(id)+1,j-lb2(id)+1,n,id) = f2d_glb(i,j,n)
              enddo
            enddo
          enddo
        enddo
      endif
!
! Do the scatter:
      len = mxsub1*mxsub2*nf ! max size of 2d subdomain
      call mpi_scatter(sndbuf,len,MPI_REAL8,rcvbuf,len,MPI_REAL8,root,
     |  MPI_COMM_WORLD,ier)
!
! Define output subdomains:
      do n=1,nf
        do j=lb2(mytid),ub2(mytid)
          do i=lb1(mytid),ub1(mytid)
            f2d_sub(i,j,n) = rcvbuf(i-lb1(mytid)+1,j-lb2(mytid)+1,n)     
          enddo
        enddo
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_scatter_f2d=time_scatter_f2d+(endtime-starttime)
      endif
      end subroutine mp_scatter_f2d
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
      real :: starttime,endtime
!
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats','Communication',ier)
      call vtbegin(100,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlats=time_bndlats+(endtime-starttime)
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
      real :: starttime,endtime
      integer :: ireq(4),istat(MPI_STATUS_SIZE,4),nfsave=0
!
#ifdef VT
!     call vtsymdef(100, 'mp_bndlats_f2d','Communication',ier)
      call vtbegin(100,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlats_f2d=time_bndlats_f2d+(endtime-starttime)
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
      real :: starttime,endtime
#ifdef VT
!     call vtsymdef(101, 'mp_bndlons','Communication',ier)
      call vtbegin(101,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlons=time_bndlons+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(101, 'mp_bndlons','Communication',ier)
      call vtend(101,ier)
#endif
      end subroutine mp_bndlons
!-----------------------------------------------------------------------
      subroutine mp_bndlons_f3d(f,id1,i0,i1,j0,j1,nf,iprint)
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
      integer :: n,ier,len,inext,iprev,isend0,isend1,nlats,j,i,k,
     |  irecv0,irecv1,lendat,jbeg,jend
      integer,save :: nfsave=0
      real,allocatable,save :: 
     |  sndbuf0(:,:,:,:), ! send buffer for i0  ,i0+1 (id1,2,mxlat+4,nf)
     |  sndbuf1(:,:,:,:), ! send buffer for i1-1,i1   (id1,2,mxlat+4,nf)
     |  rcvbuf0(:,:,:,:), ! recv buffer for i0-2,i0-1 (id1,2,mxlat+4,nf)
     |  rcvbuf1(:,:,:,:)  ! recv buffer for i1+1,i1+2 (id1,2,mxlat+4,nf)
      integer :: ireq(4),istat(MPI_STATUS_SIZE,4)
      real :: starttime,endtime
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f3d','Communication',ier)
      call vtbegin(102,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
     |write(6,"('mp_bndlons_f3d: id1=',i4,' mxlat+4=',i4,
     |  ' nf=',i4,' len=',i6,' lendat=',i6)") id1,mxlat+4,nf,len,
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
          do i=1,2
            do k=1,id1
              sndbuf0(k,i,j-jbeg+1,n) = f(k,i0+i-1,j,n)
              sndbuf1(k,i,j-jbeg+1,n) = f(k,i1-2+i,j,n)
            enddo
          enddo

          if (iprint > 0)
     |    write(6,"('mp_bndlons_f3d: n=',i3,' j0,1=',2i4,' j=',i3,
     |      ' sndbuf0(id1-1,:,j-jbeg+1,n)=',2e12.4,
     |      ' sndbuf1(id1-1,:,j-jbeg+1,n)=',2e12.4)") n,j0,j1,j,
     |      sndbuf0(id1-1,:,j-jbeg+1,n),sndbuf1(id1-1,:,j-jbeg+1,n)

        enddo ! j=jbeg,jend

        if (iprint >0) then
        write(6,"('bndlons_f3d: sndbuf0 field ',i2,' min,max=',2e12.4)")
     |    n,minval(sndbuf0(:,:,1:nlats,n)),
     |      maxval(sndbuf0(:,:,1:nlats,n))
        write(6,"('bndlons_f3d: sndbuf1 field ',i2,' min,max=',2e12.4)")
     |    n,minval(sndbuf1(:,:,1:nlats,n)),
     |      maxval(sndbuf1(:,:,1:nlats,n))
        endif
      enddo ! do n=1,nf
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlons_f3d=time_bndlons_f3d+(endtime-starttime)
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
      real :: starttime,endtime
#ifdef VT
!     call vtsymdef(102, 'mp_bndlons_f2d','Communication',ier)
      call vtbegin(102,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlons_f2d=time_bndlons_f2d+(endtime-starttime)
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
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Serial run of an MPI job (1 task):
      if (ntask==1) then
        call mk_polelat( 0,1,ixt)
        call mk_polelat(-1,2,ixt)
        call mk_polelat(nlat+1,nlat  ,ixt)
        call mk_polelat(nlat+2,nlat-1,ixt)
        return
      endif ! ntask==1
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

      if (mpi_timing) then
        endtime = mpi_wtime()
        time_polelats=time_polelats+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(103, 'mp_polelat','Communication',ier)
      call vtend(103,ier)
#endif
      end subroutine mp_polelats
!-----------------------------------------------------------------------
      subroutine mp_polelats_f3d(f,lev0,lev1,lon0,lon1,lat0,lat1,nf,
     |  polesign)
!
! Make latitudes j=-1,0 from j=2,1, and j=nlat+1,nlat+2 from nlat-1,nlat
! Longitude data are swapped across the poles, with a sign change if
!   necessary.
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,nf
      real,intent(in) :: polesign(nf)
      real,intent(inout) :: f(lev0:lev1,lon0:lon1,lat0-2:lat1+2,nf)
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
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Serial run of an MPI job (1 task):
!     if (ntask==1) then
!       call mk_polelat( 0,1,ixt)
!       call mk_polelat(-1,2,ixt)
!       call mk_polelat(nlat+1,nlat  ,ixt)
!       call mk_polelat(nlat+2,nlat-1,ixt)
!       return
!     endif ! ntask==1
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
        allocate(sndbuf(nlevp1,mxlon,2,nf,ntsend),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_polelats: error allocating sndbuf: ',
     |      ' mxlon=',i3,' ntsend=',i3)") mxlon,ntsend
      endif
      sndbuf = 0.
      lensend = nlevp1*mxlon*2*nf
!
! Allocate receive buffer.
! rcvbuf(:,:,:,:,:), ! recv for new lats on this task (k,i,2,nf,ntrecv)
!
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevp1,mxlon,2,nf,ntrecv),stat=ier)
        if (ier /= 0)
     |    write(6,"('>>> mp_polelats: error allocating rcvbuf: ',
     |      ' mxlon=',i3,' ntrecv=',i3)") mxlon,ntrecv
      endif
      rcvbuf = 0.
      lenrecv = nlevp1*mxlon*2*nf
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
          do if=1,nf
            ii = 0
            do i=isendto(1,it),isendto(2,it)
              ii = ii+1 
              sndbuf(:,ii,1:2,if,iit) = f(:,i,j0:j1,if)

!             write(6,"('mp_polelats_f3d: lon0,1=',2i4,' i=',i4,' j0,1='
!    |          ,2i4,' f(:,i,j0,if)=',2e12.4,' f(:,i,j1,if)=',2e12.4)") 
!    |          lon0,lon1,i,j0,j1,
!    |          minval(f(:,i,j0,if)),maxval(f(:,i,j0,if)),
!    |          minval(f(:,i,j1,if)),maxval(f(:,i,j1,if))

            enddo
          enddo ! if=1,nf
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
! Define f at new latitudes from recvbuf:
!
      j0 = -1 ; j1 = 0             ! Receive data for new S lats 0,-1
      if (mytidj == ntaskj-1) then ! Receive data for new N lats nlat+1,nlat+2
        j0 = nlat+1 ; j1 = nlat+2
      endif
      iit = 0
      do it=0,ntask-1
        if (irecvfm(1,it)/=-1.and.irecvfm(2,it)/=-1) then 
          iit = iit+1
          field_loop: do if=1,nf
!           if (polesign(if)==0.) cycle field_loop ! n2d,ne,o2p
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
              f(:,irecvlons(ii),j0,if)=
     |          rcvbuf(:,ii,2,if,iit)*polesign(if)
              f(:,irecvlons(ii),j1,if)=
     |          rcvbuf(:,ii,1,if,iit)*polesign(if)

!             write(6,"('mp_polelats_f3d: lon0,1=',2i4,' ii=',i4,
!    |          ' irecvlons(ii)=',i4,' j0,1=',2i3,
!    |          ' f(:,irecvlons(ii),j0,if)=',2e12.4)") lon0,lon1,
!    |          ii,irecvlons(ii),j0,j1,minval(f(:,irecvlons(ii),j0,if)),
!    |                                 maxval(f(:,irecvlons(ii),j0,if))

            enddo
          enddo field_loop
        endif
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_polelat_f3d=time_polelat_f3d+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(103, 'mp_polelat_f3d','Communication',ier)
      call vtend(103,ier)
#endif
      end subroutine mp_polelats_f3d
!-----------------------------------------------------------------------
      subroutine mp_gatherlons_f3d(f,k0,k1,i0,i1,j0,j1,nflds,name)
!
! Gather longitude data in a row of tasks to leftmost task in the row.
! On entry f(k0:k1,i0:i1,j0:j1,nflds) is defined for current task.
! On exit f(k0:k1,nlonp4,j0:j1,nflds) is defined for task with mytidi==0.
!
! Args:
      integer,intent(in) :: k0,k1,i0,i1,j0,j1,nflds
      real,intent(inout) :: f(k0:k1,nlonp4,j0:j1,nflds)
      character(len=*),intent(in) :: name
!
! Local:
      integer :: i,k,j,n,nlons,nlats,nlonrecv,nlevs,len,idest,isrc,ier,
     |  isend,irecv,itask,lonrecv0,lonrecv1,mtag
      real :: fmin,fmax
      real :: starttime,endtime
      real,allocatable,save ::
     |  sndbuf(:,:,:,:), ! send buffer (nlevs,mxlon,mxlat,nflds)
     |  rcvbuf(:,:,:,:)  ! recv buffer (nlevs,mxlon,mxlat,nflds)
#if defined(INTERCOMM) || defined(CISMAH)
      integer,parameter :: mxnf=7        ! 7 2d fields for CISM
#else
      integer,parameter :: mxnf=6        ! 6 2d fields
#endif

#ifdef VT
!     call vtsymdef(106, 'mp_gatherlons_f3d','Communication',ier)
      call vtbegin(106,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()

      if (nflds > mxnf) then
        write(6,"('>>> mp_gatherlons_f3d: nflds=',i4,' but cannot be ',
     |    'called with greater than mxnf=',i4)") nflds,mxnf
        call shutdown('mp_gatherlons_f3d')
      endif
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

      if (nlons > mxlon) then
        write(6,"('>>> mp_gatherlons_f3d: nlons=',i4,' is too big:',
     |    ' must be <= mxlon=',i4,' field=',a)") nlons,mxlon,name
        call shutdown('mp_gatherlons_f3d')
      endif
      if (nlats > mxlat) then
        write(6,"('>>> mp_gatherlons_f3d: nlats=',i4,' is too big:',
     |    ' must be <= mxlat=',i4,' field=',a)") nlats,mxlat,name
        call shutdown('mp_gatherlons_f3d')
      endif
      if (nlevs > nlevp1) then
        write(6,"('>>> mp_gatherlons_f3d: nlevs=',i4,' is too big:',
     |    ' must be <= nlevp1=',i4,' field=',a)") nlevs,nlevp1,name
        call shutdown('mp_gatherlons_f3d')
      endif
!
! Allocate at every call for debug:
!
      if (allocated(sndbuf)) deallocate(sndbuf)
      if (allocated(rcvbuf)) deallocate(rcvbuf)
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
              do k=k0,k1
              f(k,lonrecv0:lonrecv1,j,n) = rcvbuf(k,1:nlonrecv,j-j0+1,n)
              enddo
            enddo ! j=j0,j1
          enddo ! n=1,nflds
        enddo ! itask=1,ntaski-1
!
! If mytidi > 0, load send buffer, and send to task (0,mytidj):
      else ! mytidi /= 0
        idest = itask_table(0,mytidj)
        do n=1,nflds
          do j=j0,j1
            do k=k0,k1
              sndbuf(k,1:nlons,j-j0+1,n) = f(k,i0:i1,j,n)
            enddo
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gatherlons_f3d=time_gatherlons_f3d+(endtime-starttime)
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
      real :: starttime,endtime
      real,allocatable,save ::
     |  sndbuf(:,:,:,:), ! send buffer (nlevs,mxlon,nlats,nflds)
     |  rcvbuf(:,:,:,:)  ! recv buffer (nlevs,mxlon,nlats,nflds)
#ifdef VT
!     call vtsymdef(107, 'mp_scatterlons_f3d','Communication',ier)
      call vtbegin(107,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
!     if (name=='UN') write(6,"('scatterlons: i0,i1,nlons=',3i4,
!    |  ' k0,k1,nlevs=',3i4,' j0,j1,nlats=',3i4,' nflds=',i3,
!    |  ' name=',a)") i0,i1,nlons,k0,k1,nlevs,j0,j1,nlats,nflds,name

      if (nlons > mxlon) then
        write(6,"('>>> mp_scatterlons_f3d: nlons=',i4,' is too big:',
     |    ' must be <= mxlon=',i4,' field=',a)") nlons,mxlon,name
        call shutdown('mp_scatterlons_f3d')
      endif
      if (nlats > mxlat) then
        write(6,"('>>> mp_scatterlons_f3d: nlats=',i4,' is too big:',
     |    ' must be <= mxlat=',i4,' field=',a)") nlats,mxlat,name
        call shutdown('mp_scatterlons_f3d')
      endif
      if (nlevs > nlevp1) then
        write(6,"('>>> mp_scatterlons_f3d: nlevs=',i4,' is too big:',
     |    ' must be <= nlevp1=',i4,' field=',a)") nlevs,nlevp1,name
        call shutdown('mp_scatterlons_f3d')
      endif
!
! Allocate at every call for debug:
!
      if (allocated(sndbuf)) deallocate(sndbuf)
      if (allocated(rcvbuf)) deallocate(rcvbuf)
!     if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0)
     |  write(6,"('>>> mp_scatterlons_f3d: error allocating sndbuf.')")
!     endif
!     if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nlevs,mxlon,mxlat,nflds),stat=ier)
        if (ier /= 0)
     |  write(6,"('>>> mp_scatterlons_f3d: error allocating rcvbuf.')")
!     endif
!
      sndbuf = 0. ; rcvbuf = 0.
      len = nlevs*mxlon*mxlat*nflds
!
! If mytidi==0, send to other tasks in my row (mytidi>0,mytidj):
      if (mytidi == 0) then
        do itask=1,ntaski-1
          idest = itask_table(itask,mytidj)
          lonsend0 = tasks(idest)%lon0
          lonsend1 = tasks(idest)%lon1
          nlonsend = lonsend1-lonsend0+1
          mtag = idest+mytid
          do n=1,nflds
            do j=j0,j1
              do k=k0,k1
              sndbuf(k,1:nlonsend,j-j0+1,n) = f(k,lonsend0:lonsend1,j,n) 
              enddo
            enddo ! j=j0,j1
          enddo ! n=1,nflds
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
            do k=k0,k1
              f(k,i0:i1,j,n) = rcvbuf(k,1:nlons,j-j0+1,n)
            enddo
          enddo ! j=j0,j1
        enddo ! n=1,nflds
      endif
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_scatterlons_f3d=time_scatterlons_f3d+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(107, 'mp_scatterlons_f3d','Communication',ier)
      call vtend(107,ier)
#endif
      end subroutine mp_scatterlons_f3d
!-----------------------------------------------------------------------
      subroutine mp_periodic_f4d(ixt)
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
      real :: starttime,endtime
!
#ifdef VT
!     call vtsymdef(108, 'mp_periodic_f4d','Communication',ier)
      call vtbegin(108,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
!
! Use assignment and return in the special case of a serial run:
      if (ntask == 1) then     ! lon0,lon1 == 1,nlonp4
        do n=1,nf4d_hist       ! loop for fields on history at time ixt
          do j=lat0-2,lat1+2   ! latitude
            f4d(n)%data(:,lon0:lon0+1,j,ixt) =     ! e.g.: 1,2 <= 73,74
     |        f4d(n)%data(:,nlon+1:nlon+2,j,ixt)
            f4d(n)%data(:,lon1-1:lon1,j,ixt) =     ! e.g.: 75,76 <= 3,4
     |        f4d(n)%data(:,lon0+2:lon0+3,j,ixt)
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_periodic_f4d=time_periodic_f4d+(endtime-starttime)
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
      real :: starttime,endtime
!
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f3d','Communication',ier)
      call vtbegin(109,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
!
! Must always allocate these buffers because this routine is
! called with different numbers of fields, and different nlevs
! (e.g., its first called with lev0:lev1, then called from
! pdynamo with mlev0:mlev1). 
! 
      if (allocated(sndbuf)) deallocate(sndbuf)
      if (allocated(rcvbuf)) deallocate(rcvbuf)

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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_periodic_f3d=time_periodic_f3d+(endtime-starttime)
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
      real :: starttime,endtime
!
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f2d','Communication',ier)
      call vtbegin(109,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
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
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_periodic_f2d=time_periodic_f2d+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(109, 'mp_periodic_f2d','Communication',ier)
      call vtend(109,ier)
#endif
      end subroutine mp_periodic_f2d
!-----------------------------------------------------------------------
      subroutine mp_geo_halos(fmsub,lon0,lon1,lat0,lat1,nf)
!
! Exchange halo/ghost points between geographic grid subdomains for nf fields.
! Two halo points are set in both lon and lat dimensions.
! Longitude halos are done first, then latitude halos are done, including
!   longitude halos that were defined first).
!
! This routine is not actually called (mp_geo_halos_f3d is called from
! oplus), but am leaving it here in case it is needed later.
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,nf
      real,intent(inout) :: 
     |  fmsub(lon0-2:lon1+2,lat0-2:lat1+2,nf)
!
! Local:
      integer :: i,j,ifld,west,east,north,south,len,isend0,isend1,
     |  irecv0,irecv1,ier,nlats,istat(MPI_STATUS_SIZE,4),ireq(4),
     |  nlons
      real,dimension(2,lat1-lat0+1,nf) ::
     |  sndlon0,sndlon1,rcvlon0,rcvlon1
      real,dimension(2,(lon1+2)-(lon0-2)+1,nf) ::
     |  sndlat0,sndlat1,rcvlat0,rcvlat1
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Init send/recv buffers for lon halos:
      sndlon0 = 0. ; rcvlon0 = 0.
      sndlon1 = 0. ; rcvlon1 = 0.
!
! Identify east and west neightbors:
      west  = itask_table(mytidi-1,mytidj)
      east  = itask_table(mytidi+1,mytidj)
!
! Exchange lat0:lat1 (lat halos are not yet defined):
      nlats = lat1-lat0+1
      len = 2*nlats*nf

!     write(6,"('mp_geo_halos: mytid=',i3,' mytidi,j=',2i4,
!    |  ' west=',i3,' east=',i3,' lat0,1=',2i4,' nlats=',i4)") 
!    |  mytid,mytidi,mytidj,west,east,lat0,lat1,nlats
!
! Send lon0:lon0+1 to the west neighbor, and lon1-1:lon1 to the east.
!
      do ifld=1,nf
        do i=1,2
          sndlon0(i,:,ifld) = fmsub(lon0+i-1,lat0:lat1,ifld)
          sndlon1(i,:,ifld) = fmsub(lon1+i-2,lat0:lat1,ifld)
        enddo
      enddo ! ifld=1,nf
!
! Send lon0:lon0+1 to the west:
      call mpi_isend(sndlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_geo_halos send lon0:lon0+1 to west')
!
! Send lon1-1:lon1 to the east:
      call mpi_isend(sndlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |    isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_geo_halos send lon1-1:lon1 to east')
!
! Recv lon0-2:lon0-1 from west:
      call mpi_irecv(rcvlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos recv lon0-2:lon0-1 from west')
!
! Recv lon1+1:lon1+2 from east:
      call mpi_irecv(rcvlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos recv lon1+1:lon1+2 from east')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_geo_halos waitall for lons')
!
! Copy lon0-2:lon0-1 from rcvlon0, and lon1+1:lon1+2 from rcvlon1:
      do ifld=1,nf
        do i=1,2
          fmsub(lon0-3+i,lat0:lat1,ifld) = rcvlon0(i,:,ifld)  
          fmsub(lon1+i,lat0:lat1,ifld)   = rcvlon1(i,:,ifld)  
        enddo
!
! Fix special case of 2 tasks in longitude dimension:
        if (east == west) then
          do i=1,2
            fmsub(lon0-3+i,lat0:lat1,ifld) = rcvlon1(i,:,ifld)  
            fmsub(lon1+i,lat0:lat1,ifld) = rcvlon0(i,:,ifld)  
          enddo
        endif

!       if (ifld==1) then
!         do i=1,2
!           write(6,"('mp_geo_halos: i=',i2,' rcv lon0-3+i=',i4,
!    |        ' from west=',i4,' ifld=',i3,
!    |        ' fmsub(lon0-3+i,lat0:lat1,ifld)=',/,(6e12.4))")
!    |        i,lon0-3+i,west,ifld,fmsub(lon0-3+i,lat0:lat1,ifld)
!           write(6,"('mp_geo_halos: i=',i2,' rcv lon1+i=',i4,
!    |        ' from east=',i4,' ifld=',i3,
!    |        ' fmsub(lon1+i,lat0:lat1,ifld)=',/,(6e12.4))")
!    |        i,lon1+i,east,ifld,fmsub(lon1+i,lat0:lat1,ifld)
!         enddo
!       endif

      enddo ! ifld=1,nf
!
! Now exchange latitudes:
      sndlat0 = 0. ; rcvlat0 = 0.
      sndlat1 = 0. ; rcvlat1 = 0.

      south = itask_table(mytidi,mytidj-1)  ! neighbor to south
      north = itask_table(mytidi,mytidj+1)  ! neighbor to north
!
! Include halo longitudes that were defined by the exchanges above:
      nlons = (lon1+2)-(lon0-2)+1 
      len = 2*nlons*nf

!     write(6,"('mp_geo_halos: mytid=',i3,' mytidi,j=',2i4,
!    |  ' south=',i3,' north=',i3,' lon0,1=',2i4,' nlons=',i4)") 
!    |  mytid,mytidi,mytidj,south,north,lon0,lon1,nlons
!
! Send lat0:lat0+1 to south neighbor, and lat1-1:lat1 to north:
      do ifld=1,nf
        do i=1,2
          sndlat0(i,:,ifld) = fmsub(:,lat0+i-1,ifld)
          sndlat1(i,:,ifld) = fmsub(:,lat1+i-2,ifld)
        enddo
      enddo
!
! Send lat0:lat0+1 to south:
      call mpi_isend(sndlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos send lat0:lat0+1 to south')
!
! Send lat1-1:lat1 to north:
      call mpi_isend(sndlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos send lat1-1:lat1 to north')
!
! Recv lat0-2:lat0-1 from south:
      call mpi_irecv(rcvlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos recv lat0-2:lat0-1 from south')
!
! Recv lat1+1:lat1+2 from north:
      call mpi_irecv(rcvlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos recv lat1+1:lat1+2 from north')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_geo_halos waitall for lats')
!
! Copy lat0-2:lat0-1 from rcvlat0, and lat1+1:lat1+2 from rcvlat1:
      do ifld=1,nf
        do i=1,2
          fmsub(:,lat0-3+i,ifld) = rcvlat0(i,:,ifld)  
          fmsub(:,lat1+i,ifld)   = rcvlat1(i,:,ifld)  
        enddo
!
! Fix special case of 2 tasks in latitude dimension:
        if (north == south) then
          do i=1,2
            fmsub(:,lat0-3+i,ifld) = rcvlat1(i,:,ifld)  
            fmsub(:,lat1+i,ifld)   = rcvlat0(i,:,ifld)  
          enddo
        endif

!       if (ifld==1) then
!         do i=1,2
!           write(6,"('mp_geo_halos: i=',i2,' rcv lat0-3+i=',i4,
!    |        ': ifld=',i3,' fmsub(:,lat0-3+i,ifld)=',/,(6e12.4))")
!    |        i,lat0-3+i,ifld,fmsub(:,lat0-3+i,ifld)
!           write(6,"('mp_geo_halos: i=',i2,' rcv lat1+i=',i4,': ifld=',
!    |        i3,' fmsub(:,lat1+i,ifld)=',/,(6e12.4))")
!    |        i,lat1+i,ifld,fmsub(:,lat1+i,ifld)
!         enddo
!       endif

      enddo ! ifld=1,nf

      if (mpi_timing) then
        endtime = mpi_wtime()
        time_geo_halos=time_geo_halos+(endtime-starttime)
      endif
      end subroutine mp_geo_halos
!-----------------------------------------------------------------------
      subroutine mp_geo_halos_f3d(fmsub,lev0,lev1,lon0,lon1,lat0,lat1,
     |  nf)
!
! Exchange halo/ghost points between geographic grid subdomains for nf fields.
! Two halo points are set in both lon and lat dimensions.
! Longitude halos are done first, then latitude halos are done, including
!   longitude halos that were defined first).
!
! This is called from oplus to set halo points in 3d field N2.
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,nf
      real,intent(inout) :: 
     |  fmsub(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,nf)
!
! Local:
      integer :: k,i,j,ifld,west,east,north,south,len,isend0,isend1,
     |  irecv0,irecv1,ier,nlats,istat(MPI_STATUS_SIZE,4),ireq(4),
     |  nlons
      real,dimension(lev0:lev1,2,lat1-lat0+1,nf) ::
     |  sndlon0,sndlon1,rcvlon0,rcvlon1
      real,dimension(lev0:lev1,2,(lon1+2)-(lon0-2)+1,nf) ::
     |  sndlat0,sndlat1,rcvlat0,rcvlat1
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Init send/recv buffers for lon halos:
      sndlon0 = 0. ; rcvlon0 = 0.
      sndlon1 = 0. ; rcvlon1 = 0.
!
! Identify east and west neightbors:
      west  = itask_table(mytidi-1,mytidj)
      east  = itask_table(mytidi+1,mytidj)
!
! Exchange lat0:lat1 (lat halos are not yet defined):
      nlats = lat1-lat0+1
      len = (lev1-lev0+1)*2*nlats*nf

!     write(6,"('mp_geo_halos_f3d: mytid=',i3,' mytidi,j=',2i4,
!    |  ' west=',i3,' east=',i3,' lat0,1=',2i4,' nlats=',i4)") 
!    |  mytid,mytidi,mytidj,west,east,lat0,lat1,nlats
!
! Send lon0:lon0+1 to the west neighbor, and lon1-1:lon1 to the east.
!
      do ifld=1,nf
        do i=1,2
          do k=lev0,lev1
            sndlon0(k,i,:,ifld) = fmsub(k,lon0+i-1,lat0:lat1,ifld)
            sndlon1(k,i,:,ifld) = fmsub(k,lon1+i-2,lat0:lat1,ifld)
          enddo
        enddo
      enddo ! ifld=1,nf
!
! Send lon0:lon0+1 to the west:
      call mpi_isend(sndlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |  'mp_geo_halos_f3d send lon0:lon0+1 to west')
!
! Send lon1-1:lon1 to the east:
      call mpi_isend(sndlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |    isend1,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |  'mp_geo_halos_f3d send lon1-1:lon1 to east')
!
! Recv lon0-2:lon0-1 from west:
      call mpi_irecv(rcvlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d recv lon0-2:lon0-1 from west')
!
! Recv lon1+1:lon1+2 from east:
      call mpi_irecv(rcvlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d recv lon1+1:lon1+2 from east')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |  'mp_geo_halos_f3d waitall for lons')
!
! Copy lon0-2:lon0-1 from rcvlon0, and lon1+1:lon1+2 from rcvlon1:
      do ifld=1,nf
        do i=1,2
          do k=lev0,lev1
            fmsub(k,lon0-3+i,lat0:lat1,ifld) = rcvlon0(k,i,:,ifld)  
            fmsub(k,lon1+i,lat0:lat1,ifld)   = rcvlon1(k,i,:,ifld)  
          enddo
        enddo
!
! Fix special case of 2 tasks in longitude dimension:
        if (east == west) then
          do i=1,2
            do k=lev0,lev1
              fmsub(k,lon0-3+i,lat0:lat1,ifld) = rcvlon1(k,i,:,ifld)  
              fmsub(k,lon1+i,lat0:lat1,ifld) = rcvlon0(k,i,:,ifld)  
            enddo
          enddo
        endif
      enddo ! ifld=1,nf
!
! Now exchange latitudes:
      sndlat0 = 0. ; rcvlat0 = 0.
      sndlat1 = 0. ; rcvlat1 = 0.

      south = itask_table(mytidi,mytidj-1)  ! neighbor to south
      north = itask_table(mytidi,mytidj+1)  ! neighbor to north
!
! Include halo longitudes that were defined by the exchanges above:
      nlons = (lon1+2)-(lon0-2)+1 
      len = (lev1-lev0+1)*2*nlons*nf

!     write(6,"('mp_geo_halos_f3d: mytid=',i3,' mytidi,j=',2i4,
!    |  ' south=',i3,' north=',i3,' lon0,1=',2i4,' nlons=',i4)") 
!    |  mytid,mytidi,mytidj,south,north,lon0,lon1,nlons
!
! Send lat0:lat0+1 to south neighbor, and lat1-1:lat1 to north:
      do ifld=1,nf
        do i=1,2
          do k=lev0,lev1
            sndlat0(k,i,:,ifld) = fmsub(k,:,lat0+i-1,ifld)
            sndlat1(k,i,:,ifld) = fmsub(k,:,lat1+i-2,ifld)
          enddo
        enddo
      enddo
!
! Send lat0:lat0+1 to south:
      call mpi_isend(sndlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d send lat0:lat0+1 to south')
!
! Send lat1-1:lat1 to north:
      call mpi_isend(sndlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d send lat1-1:lat1 to north')
!
! Recv lat0-2:lat0-1 from south:
      call mpi_irecv(rcvlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d recv lat0-2:lat0-1 from south')
!
! Recv lat1+1:lat1+2 from north:
      call mpi_irecv(rcvlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,
     |    'mp_geo_halos_f3d recv lat1+1:lat1+2 from north')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |  'mp_geo_halos_f3d waitall for lats')
!
! Copy lat0-2:lat0-1 from rcvlat0, and lat1+1:lat1+2 from rcvlat1:
      do ifld=1,nf
        do i=1,2
          do k=lev0,lev1
            fmsub(k,:,lat0-3+i,ifld) = rcvlat0(k,i,:,ifld)  
            fmsub(k,:,lat1+i,ifld)   = rcvlat1(k,i,:,ifld)  
          enddo
        enddo
!
! Fix special case of 2 tasks in latitude dimension:
        if (north == south) then
          do i=1,2
            do k=lev0,lev1
              fmsub(k,:,lat0-3+i,ifld) = rcvlat1(k,i,:,ifld)  
              fmsub(k,:,lat1+i,ifld)   = rcvlat0(k,i,:,ifld)  
            enddo 
          enddo
        endif
      enddo ! ifld=1,nf

      if (mpi_timing) then
        endtime = mpi_wtime()
        time_geo_halos_f3d=time_geo_halos_f3d+(endtime-starttime)
      endif
      end subroutine mp_geo_halos_f3d
!-----------------------------------------------------------------------
      subroutine mp_close
      integer :: ier
      call mpi_finalize(ier)
      if (ier /= 0) then
        write(6,"(/,'>>> WARNING: error from mp_close: ier=',i3)") ier
      endif
      end subroutine mp_close
!-----------------------------------------------------------------------
      subroutine mp_mageq(fin,fout,nf,mlon0,mlon1,mlat0,mlat1,
     |  mlev0,mlev1)
!
! Each task needs values of conductivities and adotv1,2 fields at the
! at the mag equator for its longitude subdomain (and all levels), for 
! the fieldline integrations.
!
! On input, fin is ped_mag, hal_mag, adotv1_mag, adotv2_mag
!   on (i,j,k) magnetic subdomain.
! On output, fout(mlon0:mlon1,mlev0:mlev1,nf) is ped_meq, hal_meq, 
!   adotv1_meq,adotv2_meq at mag equator at longitude subdomain and 
!   all levels.
!
! Args:
      integer :: mlon0,mlon1,mlat0,mlat1,mlev0,mlev1,nf
      real,intent(in)  :: fin (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,nf)
      real,intent(out) :: fout(mlon0:mlon1,mlev0:mlev1,nf)
!
! Local:
      real,allocatable,save :: sndbuf(:,:,:), ! mxmaglon,mlev0:mlev1,nf
     |                         rcvbuf(:,:,:)  ! mxmaglon,mlev0:mlev1,nf
      integer :: i,j,n,itask,ier,len,jlateq,ireqsend,ireqrecv
      real :: starttime,endtime
      logical :: have_eq

      if (mpi_timing) starttime = mpi_wtime()
!
! Allocate buffers once per run:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxmaglon,mlev0:mlev1,nf),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_mageq: error allocating sndbuf')")
          call shutdown('mp_mageq')
        endif
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxmaglon,mlev0:mlev1,nf),stat=ier)
        if (ier /= 0) then
          write(6,"('>>> mp_mageq: error allocating rcvbuf')")
          call shutdown('mp_mageq')
        endif
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = mxmaglon*(mlev1-mlev0+1)*nf 
!
! If mag equator is in current subdomain, load it into sndbuf
! and send to other tasks in my task column (mytidi) 
!
      jlateq = (nmlat+1)/2   ! lat index of mag equator (49)
      have_eq = .false.
      do j=mlat0,mlat1
        if (j == jlateq) then ! load send buffer w/ data at equator
          have_eq = .true.
          do i=mlon0,mlon1
            sndbuf(i-mlon0+1,:,:) = fin(i,j,:,:)
          enddo
!
! Send mag equator data to other tasks in my task column (mytidi): 
          do itask=0,ntask-1
            if (itask /= mytid.and.tasks(itask)%mytidi==mytidi) then
              call mpi_isend(sndbuf,len,MPI_REAL8,itask,1,
     |          MPI_COMM_WORLD,ireqsend,ier)
              if (ier /= 0) call handle_mpi_err(ier,'mp_mageq isend')
              call mpi_wait(ireqsend,irstat,ier)
            endif ! another task in mytidi
          enddo ! itask=0,ntask-1
        endif ! j==jlateq
      enddo ! j=mlat0,mlat1
!
! Receive by other tasks in the sending task's column:
      fout = 0.
      if (.not.have_eq) then ! find task to receive from
        do itask=0,ntask-1
          do j=tasks(itask)%mlat0,tasks(itask)%mlat1
            if (j == jlateq.and.tasks(itask)%mytidi==mytidi) then
              call mpi_irecv(rcvbuf,len,MPI_REAL8,itask,1,
     |          MPI_COMM_WORLD,ireqrecv,ier)
              if (ier /= 0) call handle_mpi_err(ier,'mp_mageq irecv')
              call mpi_wait(ireqrecv,irstat,ier)
              do n=1,nf
                do i=mlon0,mlon1
                  fout(i,:,n) = rcvbuf(i-mlon0+1,:,n)

!                 write(6,"('Receive task: mxmaglon=',i4,' mlon0,1=',
!    |              2i4,' itask mlon0,1=',2i4,' n=',i4,' i=',i4,
!    |              ' fout(i,:,n)=',/,(6e12.4))") mxmaglon,mlon0,mlon1,
!    |              tasks(itask)%mlon0,tasks(itask)%mlon1,
!    |              n,i,fout(i,:,n)

                enddo
              enddo
            endif ! itask has mag eq and is in my column (sending task)
          enddo ! scan itask latitudes
        enddo ! task table search
!
! If I am the sending task, set fout to equator values of input array:
      else
        do n=1,nf
          do i=mlon0,mlon1
            fout(i,:,n) = fin(i,jlateq,:,n)

!           write(6,"('Sending task: mxmaglon=',i4,' mlon0,1=',2i4,
!    |        ' n=',i4,' i=',i4,' fout(i,:,n)=',/,(6e12.4))") mxmaglon,
!    |        mlon0,mlon1,n,i,fout(i,:,n)

          enddo
        enddo
      endif ! I am receiving or sending task
!     do n=1,nf
!       write(6,"('end mp_mageq: n=',i3,' fout min,max=',
!    |    2e12.4)") n,minval(fout(:,:,n)),maxval(fout(:,:,n))
!     enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mageq=time_mageq+(endtime-starttime)
      endif
      end subroutine mp_mageq
!-----------------------------------------------------------------------
      subroutine mp_mageq_jpm1(f,mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  feq_jpm1,nf)
!
! All tasks need data at mag latitudes equator-1, equator+1 at global 
!   longitudes. 
! On input: f is 6 fields on mag subdomains: zigm11,zigm22,zigmc,zigm2,rim1,rim2
! On output: feq_jpm1(nmlonp1,2,nf)
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nf
      real,intent(in) :: f(mlon0:mlon1,mlat0:mlat1,nf)
      real,intent(out) :: feq_jpm1(nmlonp1,2,nf) ! eq-1,eq+1
!
! Local:
      integer :: ifld,j,ier,len,jlateq,itask,i0,i1
      real,allocatable,save :: sndbuf(:,:,:,:), rcvbuf(:,:,:,:)
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Allocate buffers once per run:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nmlonp1,2,nf,0:ntask-1),stat=ier)
        if (ier /= 0) 
     |    call shutdown('mp_mageq_jpm1 error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nmlonp1,2,nf,0:ntask-1),stat=ier)
        if (ier /= 0) 
     |    call shutdown('mp_mageq_jpm1 error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = nmlonp1*2*nf
!
! Load send buffer w/ eq +/- 1 for current subdomain
! (redundant to all tasks for alltoall)
!
      jlateq = (nmlat+1)/2 
      do itask=0,ntask-1
        do j=mlat0,mlat1
          if (j == jlateq+1) then     ! equator+1
            sndbuf(mlon0:mlon1,1,:,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq-1) then ! equator-1
            sndbuf(mlon0:mlon1,2,:,itask) = f(mlon0:mlon1,j,:)
          endif ! j==jlateq
        enddo ! j=mlat0,mlat1
      enddo ! itask
!
! Do the exchange:
      call mpi_alltoall(sndbuf,len,MPI_REAL8,
     |                  rcvbuf,len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_mageq_jpm1 call mpi_alltoall')
!
! Unpack receive buffer to output feq_jpm1(nmlonp1,2,nf)   
! rcvbuf(nmlonp1,2,nf,0:ntask-1)
      do itask=0,ntask-1
        i0 = tasks(itask)%mlon0 ; i1 = tasks(itask)%mlon1
        do j=tasks(itask)%mlat0,tasks(itask)%mlat1
          if (j == jlateq+1) then     ! equator+1
            feq_jpm1(i0:i1,1,:) = rcvbuf(i0:i1,1,:,itask)
          elseif (j == jlateq-1) then ! equator-1
            feq_jpm1(i0:i1,2,:) = rcvbuf(i0:i1,2,:,itask)
          endif
        enddo ! j=mlat0,mlat1
      enddo ! itask
!
! Periodic point:
      feq_jpm1(nmlonp1,:,:) = feq_jpm1(1,:,:)

! real,intent(out) :: feq_jpm1(nmlonp1,2,nf)   
!     do ifld=1,nf
!       do j=1,2
!         write(6,"('mp_mageq_jpm1: ifld=',i3,' j=',i2,
!    |      ' feq_jpm1(:,j,ifld)=',/,(6e12.4))") ifld,
!    |      j,feq_jpm1(:,j,ifld)
!       enddo
!     enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mageq_jpm1=time_mageq_jpm1+(endtime-starttime)
      endif
      end subroutine mp_mageq_jpm1
!-----------------------------------------------------------------------
      subroutine mp_mageq_jpm3(f,mlon0,mlon1,mlat0,mlat1,nmlonp1,
     |  feq_jpm3,nf)
!
! All tasks need global longitudes at mag latitudes equator, 
!   and equator +/- 1,2,3
! On input: f is nf fields on mag subdomains
! On output: feq_jpm3(nmlonp1,-3:3,nf) has global lons at eq, eq +/- 1,2,3
!   2nd dimension of feq_jpm3 (and send/recv buffers) is as follows:
!   +3: eq+3
!   +2: eq+2
!   +1: eq+1
!    0: eq
!   -1: eq-1
!   -2: eq-2
!   -3: eq-3
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nf
      real,intent(in) :: f(mlon0:mlon1,mlat0:mlat1,nf)
      real,intent(out) :: feq_jpm3(nmlonp1,-3:3,nf)
!
! Local:
      integer :: ifld,j,ier,len,jlateq,itask,i0,i1
      real,allocatable,save :: sndbuf(:,:,:,:), rcvbuf(:,:,:,:)
      real :: starttime,endtime
#if defined(INTERCOMM) || defined(CISMAH)
      integer,parameter :: mxnf=7        ! 7 2d fields for CISM
#else
      integer,parameter :: mxnf=6        ! 6 2d fields
#endif

      if (mpi_timing) starttime = mpi_wtime()

      if (nf > mxnf) then
        write(6,"('>>> mp_mageq_jpm3: nf=',i4,' but cannot be ',
     |    'called with greater than mxnf=',i4)") nf,mxnf
        call shutdown('mp_mageq_jpm3')
      endif
!
! Allocate buffers once per run (max number of fields):
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nmlonp1,-3:3,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0) 
     |    call shutdown('mp_mageq_jpm3 error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nmlonp1,-3:3,mxnf,0:ntask-1),stat=ier) 
        if (ier /= 0) 
     |    call shutdown('mp_mageq_jpm3 error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = nmlonp1*7*nf
!
! Load send buffer w/ eq +/- 1 for current subdomain
! (redundant to all tasks for alltoall)
!
      jlateq = (nmlat+1)/2 
      do itask=0,ntask-1
        do j=mlat0,mlat1
          if (j == jlateq-3) then       ! equator-3
            sndbuf(mlon0:mlon1,-3,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq-2) then   ! equator-2
            sndbuf(mlon0:mlon1,-2,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq-1) then   ! equator-1
            sndbuf(mlon0:mlon1,-1,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq) then     ! equator
            sndbuf(mlon0:mlon1,0,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq+1) then   ! equator
            sndbuf(mlon0:mlon1,1,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq+2) then   ! equator
            sndbuf(mlon0:mlon1,2,1:nf,itask) = f(mlon0:mlon1,j,:)
          elseif (j == jlateq+3) then   ! equator
            sndbuf(mlon0:mlon1,3,1:nf,itask) = f(mlon0:mlon1,j,:)
          endif ! j==jlateq
        enddo ! j=mlat0,mlat1
      enddo ! itask
!
! Do the exchange:
      call mpi_alltoall(sndbuf(:,:,1:nf,:),len,MPI_REAL8,
     |  rcvbuf(:,:,1:nf,:),len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_mageq_jpm3 call mpi_alltoall')
!
! Unpack receive buffer to output feq_jpm3(nmlonp1,-3:3,nf)   
! rcvbuf(nmlonp1,-3:3,nf,0:ntask-1)
      do itask=0,ntask-1
        i0 = tasks(itask)%mlon0 ; i1 = tasks(itask)%mlon1
        do j=tasks(itask)%mlat0,tasks(itask)%mlat1
          if (j == jlateq-3) then         ! equator-3
            feq_jpm3(i0:i1,-3,:) = rcvbuf(i0:i1,-3,1:nf,itask)
          elseif (j == jlateq-2) then     ! equator-2
            feq_jpm3(i0:i1,-2,:) = rcvbuf(i0:i1,-2,1:nf,itask)
          elseif (j == jlateq-1) then     ! equator-1
            feq_jpm3(i0:i1,-1,:) = rcvbuf(i0:i1,-1,1:nf,itask)
          elseif (j == jlateq) then       ! equator
            feq_jpm3(i0:i1,0,:) = rcvbuf(i0:i1,0,1:nf,itask)
          elseif (j == jlateq+1) then     ! equator+1
            feq_jpm3(i0:i1,1,:) = rcvbuf(i0:i1,1,1:nf,itask)
          elseif (j == jlateq+2) then     ! equator+2
            feq_jpm3(i0:i1,2,:) = rcvbuf(i0:i1,2,1:nf,itask)
          elseif (j == jlateq+3) then     ! equator+3
            feq_jpm3(i0:i1,3,:) = rcvbuf(i0:i1,3,1:nf,itask)
          endif
        enddo ! j=mlat0,mlat1
      enddo ! itask
!
! Periodic point:
      feq_jpm3(nmlonp1,:,:) = feq_jpm3(1,:,:)

! real,intent(out) :: feq_jpm3(nmlonp1,2,nf)   
!     do ifld=1,nf
!       do j=1,2
!         write(6,"('mp_mageq_jpm3: ifld=',i3,' j=',i2,
!    |      ' feq_jpm3(:,j,ifld)=',/,(6e12.4))") ifld,
!    |      j,feq_jpm3(:,j,ifld)
!       enddo
!     enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mageq_jpm3=time_mageq_jpm3+(endtime-starttime)
      endif
      end subroutine mp_mageq_jpm3
!-----------------------------------------------------------------------
      subroutine mp_magpole_2d(f,ilon0,ilon1,ilat0,ilat1,
     |  nglblon,jspole,jnpole,fpole_jpm2,nf)
!
! Obtain global longitudes of 2d fields at highest mag latitude bands. 
! These will be used to calculate values at magnetic poles.
! Return fpole_jpm2(nglblon,1->4,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields. Buffers are allocated once per run for
! mxnf fields, then passed and processed for nf fields.
!
! Args:
      integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon,
     |  jspole,jnpole,nf
      real,intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
      real,intent(out) :: fpole_jpm2(nglblon,4,nf)
!
! Local:
      integer :: j,k,ier,len,itask,i0,i1,j0,j1
      real,allocatable,save :: sndbuf(:,:,:,:),rcvbuf(:,:,:,:)
      real :: starttime,endtime
#if defined(INTERCOMM) || defined(CISMAH)
      integer,parameter :: mxnf=7        ! 7 2d fields for CISM
#else
      integer,parameter :: mxnf=6        ! 6 2d fields
#endif

      if (mpi_timing) starttime = mpi_wtime()

      if (nf > mxnf) then
        write(6,"('>>> mp_magpole_2d: nf=',i4,' but cannot be ',
     |    'called with greater than mxnf=',i4)") nf,mxnf
        call shutdown('mp_magpole_2d')
      endif
!
! Allocate the max number of fields, only once per run:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nglblon,4,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpole_2d: error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nglblon,4,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpole_2d: error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = nglblon*4*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
! (redundant to all tasks for alltoall)
!
      do itask=0,ntask-1
        do j=ilat0,ilat1
          if (j==jspole+1) then             ! south pole +1
            sndbuf(ilon0:ilon1,1,1:nf,itask) = f(ilon0:ilon1,j,:)
          elseif (j==jspole+2) then         ! south pole +2
            sndbuf(ilon0:ilon1,2,1:nf,itask) = f(ilon0:ilon1,j,:)
          elseif (j==jnpole-1) then         ! north pole -1
            sndbuf(ilon0:ilon1,3,1:nf,itask) = f(ilon0:ilon1,j,:)
          elseif (j==jnpole-2) then         ! north pole -2
            sndbuf(ilon0:ilon1,4,1:nf,itask) = f(ilon0:ilon1,j,:)
          endif
        enddo
      enddo
!
! Do the exchange:
      call mpi_alltoall(sndbuf(:,:,1:nf,:),len,MPI_REAL8,
     |  rcvbuf(:,:,1:nf,:),len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_magpole_2d call mpi_alltoall')
!
! Unpack rcvbuf into output fpole_jpm2(nglblon,4,nf)
      do itask=0,ntask-1
        i0 = tasks(itask)%mlon0
        i1 = tasks(itask)%mlon1
        j0 = tasks(itask)%mlat0
        j1 = tasks(itask)%mlat1
        do j=j0,j1
          if (j==jspole+1) then             ! south pole +1
            fpole_jpm2(i0:i1,1,:) = rcvbuf(i0:i1,1,1:nf,itask)
          elseif (j==jspole+2) then         ! south pole +2
            fpole_jpm2(i0:i1,2,:) = rcvbuf(i0:i1,2,1:nf,itask)  
          elseif (j==jnpole-1) then         ! north pole -1
            fpole_jpm2(i0:i1,3,:) = rcvbuf(i0:i1,3,1:nf,itask)  
          elseif (j==jnpole-2) then         ! north pole -2
            fpole_jpm2(i0:i1,4,:) = rcvbuf(i0:i1,4,1:nf,itask)  
          endif
        enddo
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_magpole_2d=time_magpole_2d+(endtime-starttime)
      endif
      end subroutine mp_magpole_2d
!-----------------------------------------------------------------------
      subroutine mp_magpole_3d(f,ilon0,ilon1,ilat0,ilat1,nlev,
     |  nglblon,jspole,jnpole,fpole_jpm2,nf)
!
! Obtain global longitudes at highest mag latitude bands. 
! These will be used to calculate values at magnetic poles.
! 
! Return fpole_jpm2(nglblon,1->4,nlev,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields. Buffers are allocated once per run for
! mxnf fields, then passed and processed for nf fields.
!
! Args:
      integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon,
     |  jspole,jnpole,nf,nlev
      real,intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nlev,nf)
      real,intent(out) :: fpole_jpm2(nglblon,4,nlev,nf)
!
! Local:
      integer :: j,k,ier,len,itask,i0,i1,j0,j1
      real,allocatable,save :: sndbuf(:,:,:,:,:),rcvbuf(:,:,:,:,:)
      real :: starttime,endtime
#if defined(INTERCOMM) || defined(CISMAH)
      integer,parameter :: mxnf=7        ! 7 2d fields for CISM
#else
      integer,parameter :: mxnf=6        ! 6 2d fields
#endif

      if (mpi_timing) starttime = mpi_wtime()

      if (nf > mxnf) then
        write(6,"('>>> mp_magpole_3d: nf=',i4,' but cannot be ',
     |    'called with greater than mxnf=',i4)") nf,mxnf
        call shutdown('mp_magpole_3d')
      endif
!
! Allocate the max number of fields, only once per run:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nglblon,4,nlev,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpole_3d: error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nglblon,4,nlev,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpole_3d: error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = nglblon*4*nlev*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
! (redundant to all tasks for alltoall)
!
      do itask=0,ntask-1
        do j=ilat0,ilat1
          do k=1,nlev
          if (j==jspole+1) then             ! south pole +1
            sndbuf(ilon0:ilon1,1,k,1:nf,itask) = f(ilon0:ilon1,j,k,:)
          elseif (j==jspole+2) then         ! south pole +2
            sndbuf(ilon0:ilon1,2,k,1:nf,itask) = f(ilon0:ilon1,j,k,:)
          elseif (j==jnpole-1) then         ! north pole -1
            sndbuf(ilon0:ilon1,3,k,1:nf,itask) = f(ilon0:ilon1,j,k,:)
          elseif (j==jnpole-2) then         ! north pole -2
            sndbuf(ilon0:ilon1,4,k,1:nf,itask) = f(ilon0:ilon1,j,k,:)
          endif
          enddo
        enddo
      enddo
!
! Do the exchange:
      call mpi_alltoall(sndbuf(:,:,:,1:nf,:),len,MPI_REAL8,
     |  rcvbuf(:,:,:,1:nf,:),len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_magpole_3d call mpi_alltoall')
!
! Unpack rcvbuf into output fpole_jpm2(nglblon,4,nf)
      do itask=0,ntask-1
        i0 = tasks(itask)%mlon0
        i1 = tasks(itask)%mlon1
        j0 = tasks(itask)%mlat0
        j1 = tasks(itask)%mlat1
        do j=j0,j1
          if (j==jspole+1) then             ! south pole +1
            fpole_jpm2(i0:i1,1,:,:) = rcvbuf(i0:i1,1,:,1:nf,itask)
          elseif (j==jspole+2) then         ! south pole +2
            fpole_jpm2(i0:i1,2,:,:) = rcvbuf(i0:i1,2,:,1:nf,itask)  
          elseif (j==jnpole-1) then         ! north pole -1
            fpole_jpm2(i0:i1,3,:,:) = rcvbuf(i0:i1,3,:,1:nf,itask)  
          elseif (j==jnpole-2) then         ! north pole -2
            fpole_jpm2(i0:i1,4,:,:) = rcvbuf(i0:i1,4,:,1:nf,itask)  
          endif
        enddo
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_magpole_3d=time_magpole_3d+(endtime-starttime)
      endif
      end subroutine mp_magpole_3d
!-----------------------------------------------------------------------
      subroutine mp_magpoles(f,ilon0,ilon1,ilat0,ilat1,
     |  nglblon,jspole,jnpole,fpoles,nf)
!
! Similiar to mp_magpole_2d, but returns global longitudes for
! j==1 and j==nmlat (not for poles +/- 2)
! Return fpoles(nglblon,2,nf) as:
!   1: j = jspole (spole)
!   2: j = jnpole (npole)
! This can be called with different number of fields nf, but cannot
! be called w/ > mxnf fields. Buffers are allocated once per run for
! mxnf fields, then passed and processed for nf fields.
!
! Args:
      integer,intent(in) :: ilon0,ilon1,ilat0,ilat1,nglblon,
     |  jspole,jnpole,nf
      real,intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
      real,intent(out) :: fpoles(nglblon,2,nf)
!
! Local:
      integer :: j,k,ier,len,itask,i0,i1,j0,j1
      real,allocatable,save :: sndbuf(:,:,:,:),rcvbuf(:,:,:,:)
      real :: starttime,endtime
      integer,parameter :: mxnf=nmlevp1

      if (mpi_timing) starttime = mpi_wtime()

      if (nf > mxnf) then
        write(6,"('>>> mp_magpoles: nf=',i4,' but cannot be ',
     |    'called with greater than mxnf=',i4)") nf,mxnf
        call shutdown('mp_magpoles')
      endif
!
! Allocate the max number of fields, only once per run:
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(nglblon,2,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpoles: error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(nglblon,2,mxnf,0:ntask-1),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_magpoles: error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      len = nglblon*2*nf
!
! Load send buffer with values at poles +/- 2 for current subdomain
! (redundant to all tasks for alltoall)
!
      do itask=0,ntask-1
        do j=ilat0,ilat1
          if (j==jspole) then             ! south pole
            sndbuf(ilon0:ilon1,1,1:nf,itask) = f(ilon0:ilon1,j,:)
          elseif (j==jnpole) then         ! npole pole
            sndbuf(ilon0:ilon1,2,1:nf,itask) = f(ilon0:ilon1,j,:)
          endif
        enddo
      enddo
!
! Do the exchange:
      call mpi_alltoall(sndbuf(:,:,1:nf,:),len,MPI_REAL8,
     |  rcvbuf(:,:,1:nf,:),len,MPI_REAL8,MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_magpoles call mpi_alltoall')
!
! Unpack rcvbuf into output fpoles(nglblon,2,nf)
      do itask=0,ntask-1
        i0 = tasks(itask)%mlon0
        i1 = tasks(itask)%mlon1
        j0 = tasks(itask)%mlat0
        j1 = tasks(itask)%mlat1
        do j=j0,j1
          if (j==jspole) then             ! south pole
            fpoles(i0:i1,1,:) = rcvbuf(i0:i1,1,1:nf,itask)
          elseif (j==jnpole) then         ! north pole
            fpoles(i0:i1,2,:) = rcvbuf(i0:i1,2,1:nf,itask)  
          endif
        enddo
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_magpoles=time_magpoles+(endtime-starttime)
      endif
      end subroutine mp_magpoles
!-----------------------------------------------------------------------
      subroutine conjugate_points
!
! Local:
      integer :: ier,i,j,js,jn,itask,ii,jj
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! nsend_south(ntask): number of lats in south to send north
! nrecv_north(ntask): number of lats in north to recv from south
!
      allocate(nsend_south(0:ntask-1),stat=ier)
      allocate(nrecv_north(0:ntask-1),stat=ier)
!
! send_south_coords: south j lats to send north
! recv_north_coords: north j lats to recv from south
!
      allocate(send_south_coords(mxmaglat,0:ntask-1),stat=ier)
      allocate(recv_north_coords(mxmaglat,0:ntask-1),stat=ier)

      nsend_south(:) = 0
      nrecv_north(:) = 0
      send_south_coords(:,:) = 0
      recv_north_coords(:,:) = 0

      magloop: do j=mlat0,mlat1
!
! In north hem: find tasks w/ conjugate points in south to recv:
! (nmlath is in params module)
        if (gmlat(j) > 0.) then  ! in north hem of current task
          js = nmlath-(j-nmlath) ! j index to south conjugate point (should be -j)
          do itask=0,ntask-1
            do jj = tasks(itask)%mlat0,tasks(itask)%mlat1
!
! Receive these north coords from the south:
              if (jj==js.and.mlon0==tasks(itask)%mlon0.and.
     |                       mlon1==tasks(itask)%mlon1) then
                nrecv_north(itask) = nrecv_north(itask)+1
                recv_north_coords(nrecv_north(itask),itask) = j
              endif
            enddo ! jj of remote task
          enddo ! itask=0,ntask-1
          if (all(nrecv_north==0)) 
     |    write(6,"('>>> WARNING: could not find north conjugate points'
     |      ,' corresponding to south latitude js=',i4,' gmlat(js)=',
     |      f8.2)") js,gmlat(js)
!
! In south hem: find tasks w/ conjugate points in north to send:
        elseif (gmlat(j) < 0..and.j /= nmlath) then ! in south hem
          jn = nmlath+(nmlath-j) ! j index of north conjugate point
          do itask=0,ntask-1
            do jj = tasks(itask)%mlat0,tasks(itask)%mlat1
              if (jj==jn.and.mlon0==tasks(itask)%mlon0.and.
     |                       mlon1==tasks(itask)%mlon1) then
                nsend_south(itask) = nsend_south(itask)+1
! Send these south coords to the north:
                send_south_coords(nsend_south(itask),itask) = j
              endif
            enddo ! jj of remote task
          enddo ! itask=0,ntask-1
          if (all(nsend_south==0)) 
     |    write(6,"('>>> WARNING: could not find south conjugate points'
     |      ,' corresponding to north latitude jn=',i4,' gmlat(jn)=',
     |      f8.2)") jn,gmlat(jn)
        endif ! in north or south hem
      enddo magloop ! j=mlat0,mlat1

!     do itask=0,ntask-1
!       if (nrecv_north(itask) > 0) then
!         write(6,"('I will receive ',i4,' north conjugate lats ',
!    |      'from task ',i4,':',/,'recv_north_coords =',/,(8i5))") 
!    |      nrecv_north(itask),itask,
!    |      recv_north_coords(1:nrecv_north(itask),itask)
!       endif
!       if (nsend_south(itask) > 0) then
!         write(6,"('I will send ',i4,' south conjugate lats to ',
!    |      'task ',i4,':',/,'send_south_coords =',
!    |      /,(8i5))") nsend_south(itask),itask,
!    |      send_south_coords(1:nsend_south(itask),itask)
!       endif
!     enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_conjugate_points=time_conjugate_points+(endtime-starttime)
      endif
      end subroutine conjugate_points
!-----------------------------------------------------------------------
      subroutine mp_mag_foldhem(f,mlon0,mlon1,mlat0,mlat1,nf)
!
! For each point in northern hemisphere (if any) of the current task
!   subdomain, receive data from conjugate point in the south (from the 
!   south task that owns it), and sum it to the north point data. 
!   Do this for nf fields. Conjugate point indices to send/recv to/from 
!   each task were determined by sub conjugate_points (this module).
! nsend_south,       ! number of south lats to send to north (each task)
! nrecv_north        ! number of north lats to send to south (each task)
!
! This routine is called from pdynamo at every timestep.
! Sub conjugate_points is called once per run, from mp_distribute.
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf
      real,intent(inout) :: f(mlon0:mlon1,mlat0:mlat1,nf)
!
! Local:
      integer :: i,j,n,len,itask,ifld,ier,nmlons
      real :: sndbuf(mxmaglon,mxmaglat,nf,0:ntask-1)
      real :: rcvbuf(mxmaglon,mxmaglat,nf,0:ntask-1)
      integer :: jsend(0:ntask-1),jrecv(0:ntask-1)
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
      sndbuf = 0. ; rcvbuf = 0.
      jsend  = 0  ; jrecv  = 0
      len = mxmaglon*mxmaglat*nf
      nmlons = mlon1-mlon0+1
!
! Send south data to north itask:
! (To avoid deadlock, do not send if north task is also myself. This will
!  happen when there is an odd number of tasks in the latitude dimension,
!  e.g., ntask == 12, 30, etc)
!
      do itask=0,ntask-1
        if (nsend_south(itask) > 0 .and. itask /= mytid) then
          do ifld = 1,nf
            do n=1,nsend_south(itask)
              sndbuf(1:nmlons,n,ifld,itask) =
     |          f(:,send_south_coords(n,itask),ifld)
            enddo
          enddo ! ifld=1,nf
          call mpi_isend(sndbuf(1,1,1,itask),len,MPI_REAL8,
     |      itask,1,MPI_COMM_WORLD,jsend(itask),ier)
          call mpi_wait(jsend(itask),irstat,ier)
        endif ! nsend_south(itask) > 0
      enddo ! itask=0,ntask-1
!
! Receive north data from south itask and add to north,
! i.e., north = north+south. (do not receive if south task is
! also myself, but do add south data to my north points, see below)
!
      do itask=0,ntask-1
        if (nrecv_north(itask) > 0 .and. itask /= mytid) then
          call mpi_irecv(rcvbuf(1,1,1,itask),len,MPI_REAL8,
     |      itask,1,MPI_COMM_WORLD,jrecv(itask),ier)
          call mpi_wait(jrecv(itask),irstat,ier)
          do ifld=1,nf
            do n=1,nrecv_north(itask)
!
! Receive lats in reverse order:
              f(mlon0:mlon1,
     |          recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) = 
     |        f(mlon0:mlon1,
     |          recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) +
     |        rcvbuf(1:nmlons,n,ifld,itask)
            enddo ! n=1,nrecv_north(itask)
          enddo ! ifld=1,nf
!
! If I am send *and* receive task, simply add my south data to my north points:
        elseif (nrecv_north(itask) > 0 .and. itask == mytid) then
          do ifld=1,nf
            do n=1,nrecv_north(itask)
              f(mlon0:mlon1,
     |          recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) = 
     |        f(mlon0:mlon1,
     |          recv_north_coords(nrecv_north(itask)-n+1,itask),ifld) +
     |        f(mlon0:mlon1,send_south_coords(n,itask),ifld)
            enddo ! n=1,nrecv_north(itask)
          enddo ! ifld=1,nf
        endif ! nrecv_north(itask) > 0
      enddo ! itask=0,ntask-1
!
! Mag equator is also "folded", but not included in conjugate points,
! so double it here:
      do j=mlat0,mlat1
        if (j==nmlath) then
          do ifld=1,nf
            f(:,j,ifld) = f(:,j,ifld)+f(:,j,ifld)
          enddo
        endif
      enddo
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mag_foldhem=time_mag_foldhem+(endtime-starttime)
      endif
      end subroutine mp_mag_foldhem
!-----------------------------------------------------------------------
      subroutine mp_mag_periodic_f2d(f,mlon0,mlon1,mlat0,mlat1,nf) 
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf
      real,intent(inout) :: f(mlon0:mlon1,mlat0:mlat1,nf)
!
! Local:
      integer :: j,ier,idest,isrc,len,ireqsend,ireqrecv,ifld,msgtag
      real :: sndbuf(mxmaglat,nf),rcvbuf(mxmaglat,nf)
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
      len = mxmaglat*nf
!
! I am a western-most task. Send lon 1 to eastern-most tasks:
      if (mytidi==0) then
        idest = itask_table_mag(ntaski-1,mytidj)
        do j=mlat0,mlat1
          sndbuf(j-mlat0+1,:) = f(1,j,:)
        enddo
        msgtag = mytid
        call mpi_isend(sndbuf,len,MPI_REAL8,idest,msgtag,MPI_COMM_WORLD,
     |    ireqsend,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_mag_periodic_f2d send to idest')
        call mpi_wait(ireqsend,irstat,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_mag_periodic_f2d wait for send')
!
! I am eastern-most task. Receive lon 1 from western-most tasks,
! and assign to nmlonp1:
      elseif (mytidi==ntaski-1) then
        isrc = itask_table_mag(0,mytidj)
        msgtag = isrc
        call mpi_irecv(rcvbuf,len,MPI_REAL8,isrc,msgtag,MPI_COMM_WORLD,
     |    ireqrecv,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_periodic_f2d recv fm isrc')
        call mpi_wait(ireqrecv,irstat,ier)
        if (ier /= 0)
     |    call handle_mpi_err(ier,'mp_mag_periodic_f2d wait for recv')

        do j=mlat0,mlat1
          f(nmlonp1,j,:) = rcvbuf(j-mlat0+1,:)
        enddo
      endif ! mytidi == 0 or ntaski-1
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mag_periodic_f2d=time_mag_periodic_f2d+(endtime-starttime)
      endif
      end subroutine mp_mag_periodic_f2d
!-----------------------------------------------------------------------
      subroutine mp_gather_pdyn(fmsub,mlon0,mlon1,mlat0,mlat1,
     |  fmglb,nmlonp1,nmlat,nf)
!
! Gather fields on mag subdomains to root task, so root task can
! complete non-parallel portion of dynamo (starting after rhspde)
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlonp1,nmlat,nf
      real,intent(in) :: fmsub(mlon0:mlon1,mlat0:mlat1,nf)
      real,intent(out) :: fmglb(nmlonp1,nmlat,nf)
!
! Local:
      integer :: len,i,j,ifld,ier,itask,i0,i1,j0,j1,nlons,nlats
      real,dimension(mxmaglon,mxmaglat,nf) :: sndbuf
      real,dimension(mxmaglon,mxmaglat,nf,0:ntask-1) :: rcvbuf
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
      sndbuf = 0.
      rcvbuf = 0.
      len = mxmaglon*mxmaglat*nf
      nlons = mlon1-mlon0+1
      nlats = mlat1-mlat0+1
!
! Load send buffer with my subdomain:
      do ifld=1,nf
        do j=mlat0,mlat1
          sndbuf(1:nlons,j-mlat0+1,ifld) = fmsub(mlon0:mlon1,j,ifld)
        enddo
      enddo
!
! Gather to root:
      call mpi_gather(sndbuf,len,MPI_REAL8,rcvbuf,len,MPI_REAL8,0,
     |  MPI_COMM_WORLD,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_gather_pdyn: mpi_gather to root')
!
! Root task unpack from receive buffer:
      if (mytid==0) then
        do itask=0,ntask-1
          i0 = tasks(itask)%mlon0 ; i1 = tasks(itask)%mlon1
          j0 = tasks(itask)%mlat0 ; j1 = tasks(itask)%mlat1
          do ifld=1,nf
            do j=j0,j1
              fmglb(i0:i1,j,ifld) = rcvbuf(1:i1-i0+1,j-j0+1,ifld,itask) 
            enddo ! j=mlat0,mlat1
          enddo ! ifld=1,nf
        enddo ! itask=0,ntask-1
      endif ! root task
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_gather_pdyn=time_gather_pdyn+(endtime-starttime)
      endif
      end subroutine mp_gather_pdyn
!-----------------------------------------------------------------------
      subroutine mp_allgather_2d(fsub,fglb,nlon_glb,nlat_glb,
     |  lon0,lon1,lat0,lat1,mag)
!
! Called from current.F90:
! call mp_allgather_2d(dkqphi,dkqphi_glb,nmlonp1,nmlat, &
!   mlon0,mlon1,mlat0,mlat1)
!
! Args:
      integer,intent(in) :: nlon_glb,nlat_glb,lon0,lon1,lat0,lat1,mag
      real,intent(in)  :: fsub(lon0:lon1,lat0:lat1)
      real,intent(out) :: fglb(nlon_glb,nlat_glb)
!
! Local:
      integer :: len,ier,i,j,i0,i1,j0,j1,itask
      real :: sndbuf(nlon_glb,nlat_glb)
      real :: rcvbuf(nlon_glb,nlat_glb,0:ntask-1)

! MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR)
! <type> SENDBUF(*), RECVBUF(*)
! INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, IERROR 

      sndbuf = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          sndbuf(i,j) = fsub(i,j)
        enddo
      enddo

      len = nlon_glb*nlat_glb
      call mpi_allgather(sndbuf,len,MPI_REAL8,
     |                   rcvbuf,len,MPI_REAL8,MPI_COMM_WORLD,ier) 
      if (ier /= 0) then
        call handle_mpi_err(ier,
     |    '>>> mp_allgather_2d: error return from mpi_allgather')
      else
!       write(6,"('mp_allgather_2d success: rcvbuf min,max=',2e12.4)")
!    |    minval(rcvbuf),maxval(rcvbuf)
      endif

      do itask=0,ntask-1
        if (mag > 0) then
          i0 = tasks(itask)%mlon0
          i1 = tasks(itask)%mlon1
          j0 = tasks(itask)%mlat0
          j1 = tasks(itask)%mlat1
        else
          i0 = tasks(itask)%lon0
          i1 = tasks(itask)%lon1
          j0 = tasks(itask)%lat0
          j1 = tasks(itask)%lat1
        endif
        do j=j0,j1
          do i=i0,i1
            fglb(i,j) = rcvbuf(i,j,itask)
          enddo
        enddo
      enddo

      end subroutine mp_allgather_2d
!-----------------------------------------------------------------------
      subroutine mp_mag_halos(fmsub,mlon0,mlon1,mlat0,mlat1,nf0,nf1)
!
! Exchange halo/ghost points between magnetic grid subdomains for nf fields.
! Only a single halo point is required in both lon and lat dimensions.
! Note that all tasks in any row of the task matrix have the same 
!   mlat0,mlat1, and that all tasks in any column of the task matrix
!   have the same mlon0,mlon1.
! Longitude halos are done first, exchanging mlat0:mlat1, then latitude
!   halos are done, exchanging mlon0-1:mlon1+1 (i.e., including the
!   longitude halos that were defined first).
!
! Args:
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nf0,nf1
      real,intent(inout) :: fmsub(mlon0-1:mlon1+1,mlat0-1:mlat1+1,
     |                            nf0:nf1)
!
! Local:
      integer :: i,j,ifld,west,east,north,south,len,isend0,isend1,
     |  irecv0,irecv1,ier,nmlats,istat(MPI_STATUS_SIZE,4),ireq(4),
     |  nmlons,nf
      real,dimension(mlat1-mlat0+1,nf0:nf1) ::
     |  sndlon0,sndlon1,rcvlon0,rcvlon1
      real,dimension((mlon1+1)-(mlon0-1)+1,nf0:nf1) ::
     |  sndlat0,sndlat1,rcvlat0,rcvlat1
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Init send/recv buffers for lon halos:
      sndlon0 = 0. ; rcvlon0 = 0.
      sndlon1 = 0. ; rcvlon1 = 0.
!
! Identify east and west neightbors:
! 27:forrtl: severe (408): fort: (2): Subscript #2 of the array ITASK_TABLE_MAG has value 6 which is greater than the upper bound of 5
      west  = itask_table_mag(mytidi-1,mytidj)
      east  = itask_table_mag(mytidi+1,mytidj)
!
! Exchange mlat0:mlat1 (lat halos are not yet defined):
      nf = nf1-nf0+1
      nmlats = mlat1-mlat0+1
      len = nmlats*nf

!     write(6,"('mp_mag_halos: mytid=',i3,' mytidi,j=',2i4,
!    |  ' west=',i3,' east=',i3,' mlat0,1=',2i4,' nmlats=',i4)") 
!    |  mytid,mytidi,mytidj,west,east,mlat0,mlat1,nmlats
!
! Send mlon0 to the west neighbor, and mlon1 to the east.
! However, tasks are periodic in longitude (see itask_table_mag),
!   and far west tasks send mlon0+1, and far east tasks send mlon1-1
!
      do ifld=nf0,nf1
! Far west tasks send mlon0+1 to far east (periodic) tasks:
        if (mytidi==0) then
          sndlon0(:,ifld) = fmsub(mlon0+1,mlat0:mlat1,ifld)
! Interior tasks send mlon0 to west neighbor:
        else
          sndlon0(:,ifld) = fmsub(mlon0,mlat0:mlat1,ifld)
        endif

! Far east tasks send mlon1-1 to far west (periodic) tasks:
        if (mytidi==nmagtaski-1) then 
          sndlon1(:,ifld) = fmsub(mlon1-1,mlat0:mlat1,ifld)
! Interior tasks send mlon1 to east neighbor:
        else
          sndlon1(:,ifld) = fmsub(mlon1,mlat0:mlat1,ifld)
        endif
      enddo ! ifld=nf0,nf1
!
! Send mlon0 to the west:
      call mpi_isend(sndlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos send mlon0 to west')
!
! Send mlon1 to the east:
      call mpi_isend(sndlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |    isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos send mlon1 to east')
!
! Recv mlon0-1 from west:
      call mpi_irecv(rcvlon0,len,MPI_REAL8,west,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos recv mlon0 from west')
!
! Recv mlon1+1 from east:
      call mpi_irecv(rcvlon1,len,MPI_REAL8,east,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos recv mlon1 from east')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos waitall for lons')
!
! Copy mlon0-1 from rcvlon0, and mlon1+1 from rcvlon1:
      do ifld=nf0,nf1
        fmsub(mlon0-1,mlat0:mlat1,ifld) = rcvlon0(:,ifld)  
        fmsub(mlon1+1,mlat0:mlat1,ifld) = rcvlon1(:,ifld)  
!
! Fix special case of 2 tasks in longitude dimension:
        if (east == west) then
          fmsub(mlon0-1,mlat0:mlat1,ifld) = rcvlon1(:,ifld)  
          fmsub(mlon1+1,mlat0:mlat1,ifld) = rcvlon0(:,ifld)  
        endif

!       if (ifld==1) then
!         write(6,"('mp_mag_halos: rcv mlon0-1=',i4,' from west=',i4,
!    |      ' ifld=',i3,' fmsub(mlon0-1,mlat0:mlat1,ifld)=',/,(6e12.4))
!    |      ") mlon0-1,west,ifld,fmsub(mlon0-1,mlat0:mlat1,ifld)
!         write(6,"('mp_mag_halos: rcv mlon1+1=',i4,' from east=',i4,
!    |      ' ifld=',i3,' fmsub(mlon1+1,mlat0:mlat1,ifld)=',/,(6e12.4))
!    |      ") mlon1+1,east,ifld,fmsub(mlon1+1,mlat0:mlat1,ifld)
!       endif

      enddo ! ifld=nf0,nf1
!
! Now exchange latitudes:
      sndlat0 = 0. ; rcvlat0 = 0.
      sndlat1 = 0. ; rcvlat1 = 0.

      south = itask_table_mag(mytidi,mytidj-1)  ! neighbor to south

! With 8 procs: 4-7:forrtl: severe (408): fort: (2): Subscript #2 of the array ITASK_TABLE_MAG has value 3 which is greater than the upper bound of 2

      north = itask_table_mag(mytidi,mytidj+1)  ! neighbor to north
!
! Include halo longitudes that were defined by the exchanges above:
      nmlons = (mlon1+1)-(mlon0-1)+1 
      len = nmlons*nf

!     write(6,"('mp_mag_halos: mytid=',i3,' mytidi,j=',2i4,
!    |  ' south=',i3,' north=',i3,' mlon0,1=',2i4,' nmlons=',i4)") 
!    |  mytid,mytidi,mytidj,south,north,mlon0,mlon1,nmlons
!
! Send mlat0 to south neighbor, and mlat1 to north:
      do ifld=nf0,nf1
        sndlat0(:,ifld) = fmsub(:,mlat0,ifld)
        sndlat1(:,ifld) = fmsub(:,mlat1,ifld)
      enddo
!
! Send mlat0 to south:
      call mpi_isend(sndlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  isend0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos send mlat0 to south')
!
! Send mlat1 to north:
      call mpi_isend(sndlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  isend1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos send mlat1 to north')
!
! Recv mlat0-1 from south:
      call mpi_irecv(rcvlat0,len,MPI_REAL8,south,1,MPI_COMM_WORLD,
     |  irecv0,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos recv mlat0-1 from south')
!
! Recv mlat1+1 from north:
      call mpi_irecv(rcvlat1,len,MPI_REAL8,north,1,MPI_COMM_WORLD,
     |  irecv1,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos recv mlat1+1 from north')
!
! Wait for completions:
      ireq = (/isend0,isend1,irecv0,irecv1/)
      istat = 0
      call mpi_waitall(4,ireq,istat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_mag_halos waitall for lats')
!
! Copy mlat0-1 from rcvlat0, and mlat1+1 from rcvlat1:
      do ifld=nf0,nf1
        fmsub(:,mlat0-1,ifld) = rcvlat0(:,ifld)  
        fmsub(:,mlat1+1,ifld) = rcvlat1(:,ifld)  

!       if (ifld==1) then
!         write(6,"('mp_mag_halos: rcv mlat0-1=',i4,': ifld=',i3,
!    |      ' fmsub(:,mlat0-1,ifld)=',/,(6e12.4))")
!    |      mlat0-1,ifld,fmsub(:,mlat0-1,ifld)
!         write(6,"('mp_mag_halos: rcv mlat1+1=',i4,': ifld=',i3,
!    |      ' fmsub(:,mlat1+1,ifld)=',/,(6e12.4))")
!    |      mlat1+1,ifld,fmsub(:,mlat1+1,ifld)
!       endif

      enddo ! ifld=nf0,nf1
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mag_halos=time_mag_halos+(endtime-starttime)
      endif
      end subroutine mp_mag_halos
!-----------------------------------------------------------------------
      subroutine mp_scatter_coeffs(rhs,cofum,coef_rhs,coef_cofum)
!
! Dynamo solver coefficients have been calculated globally by the
! root task in pdynamo. Now scatter these out to non-root tasks.
! The root task broadcasts the global arrays, then the non-root
! tasks sort out their own subdomains. This is called by the root
! task from pdynamo, and called by non-root tasks from pdynamo.
! (non-root tasks will be waiting here until the root task makes the call)
! Note: root task should not reference coef_rhs or coef_cofum
!
! Args:
      real,intent(in) :: rhs(nmlonp1,nmlat),cofum(nmlonp1,nmlat,9)
      real,intent(out) :: coef_rhs(mlon0:mlon1,mlat0:mlat1)
      real,intent(out) :: coef_cofum(mlon0:mlon1,mlat0:mlat1,9)
!
! Local:
      integer :: ier,len,i,j,n
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Broadcast rhs (from pdynamo c0(:,:,10))
      len = nmlat*nmlonp1
      call mpi_bcast(rhs,len,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |    'mp_scatter_coeffs: bcast global rhs')
!     write(6,"('mp_scatter_coeffs after bcast: rhs min,max=',2e12.4)")
!    |  minval(rhs),maxval(rhs)
!
! Broadcast cofum (from pdynamo cofum(:,:,1-9)):
      len = nmlat*nmlonp1*9
      call mpi_bcast(cofum,len,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |    'mp_scatter_coeffs: bcast global cofum')
!     write(6,"('mp_scatter_coeffs after bcast: cofum min,max=',
!    |  2e12.4)") minval(cofum),maxval(cofum)

      do j=mlat0,mlat1
        do i=mlon0,mlon1
          coef_rhs(i,j) = rhs(i,j)
          coef_cofum(i,j,:) = cofum(i,j,:)
        enddo
      enddo
!     write(6,"('mp_scatter_coeffs: coef_rhs   min,max=',2e12.4)")
!    |  minval(coef_rhs),maxval(coef_rhs)
!     write(6,"('mp_scatter_coeffs: coef_cofum min,max=',2e12.4)")
!    |  minval(coef_cofum),maxval(coef_cofum)

      if (mpi_timing) then
        endtime = mpi_wtime()
        time_scatter_coeffs=time_scatter_coeffs+(endtime-starttime)
      endif
      end subroutine mp_scatter_coeffs
!-----------------------------------------------------------------------
      subroutine mp_scatter_phim(phim_glb,phim)
      real,intent(in) :: phim_glb(nmlonp1,nmlat)
      real,intent(out) :: phim(mlon0:mlon1,mlat0:mlat1)
!
! Local:
      integer :: ier,len,i,j
      real :: starttime,endtime

      if (mpi_timing) starttime = mpi_wtime()
!
! Broadcast global phim (from pdynamo phim(nmlonp1,nmlat)):
      len = nmlat*nmlonp1
      call mpi_bcast(phim_glb,len,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      if (ier /= 0) call handle_mpi_err(ier,
     |    'mp_scatter_phim: bcast global phim')
!     write(6,"('mp_scatter_phim after bcast: phim min,max=',2e12.4)")
!    |  minval(phim),maxval(phim)
!
! Define subdomains:
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          phim(i,j) = phim_glb(i,j)
        enddo
      enddo
!     write(6,"('mp_scatter_phim: phim min,max=',2e12.4)")
!    |  minval(phim),maxval(phim)
      
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_scatter_phim=time_scatter_phim+(endtime-starttime)
      endif
      end subroutine mp_scatter_phim
!-----------------------------------------------------------------------
      subroutine mp_bcast_real(buff,len,root)
!
! Args:
      integer,intent(in) :: len,root
      real,intent(inout) :: buff(len)
!
! Local:
      integer :: ier
!
!int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, 
!               MPI_Comm comm )

      call mpi_bcast(buff,len,MPI_REAL8,root,MPI_COMM_WORLD,ier)
      if (ier /= 0) then
        call handle_mpi_err(ier,
     |    '>>> mp_bcast_real: error return from mpi_bcast')
      else
!       write(6,"('mp_bcast_real success: len=',i12,' root=',i8,
!    |    ' buff=',(6e12.4))") len,root,buff
      endif
!
      end subroutine mp_bcast_real
!-----------------------------------------------------------------------
      subroutine mp_mag_jslot(fin,mlon00,mlon11,mlat00,mlat11,
     |  fout,jneed,mxneed,nf)
!
! Current task needs to receive (from other tasks) field f at (non-zero) 
! latitude indices in jneed, at all longitudes in the current subdomain.
! Note subdomains include halo points mlon0-1 and mlat1+1. Data in f also
! includes halo points (will need the lat data at halo-longitudes)
!
! Args:
      integer,intent(in) :: mlon00,mlon11,mlat00,mlat11 ! subdomains w/ halos
      integer,intent(in) :: nf ! number of fields
      integer,intent(in) ::  mxneed ! max number of needed lats (nmlat+2)
      integer,intent(in) :: jneed(mxneed) ! j-indices of needed lats (where /= -1)
      real,intent(in) :: fin(mlon00:mlon11,mlat00:mlat11,nf) ! data at current subdomain
      real,intent(out) :: fout(mlon00:mlon11,mxneed,nf) ! returned data at needed lats
!
! Local:
      integer :: ier,njneed,i,j,n,nj,idest,ireq(2),
     |  isend(ntask),irecv(ntask),icount,
     |  len,istat,nlons,isrc,msgid,ifld
      integer :: tij ! rank in cols_comm (0 to nmagtaskj-1)
      integer :: jhave(mxneed),njhave,wid
      integer :: peersneed(mxneed,0:nmagtaskj-1)
      integer :: jneedall (mxneed,0:nmagtaskj-1)
      real,allocatable :: sndbuf(:,:,:)
      real,allocatable :: rcvbuf(:,:,:)
      real :: starttime,endtime
      integer,external :: ixfind

      if (mpi_timing) starttime = mpi_wtime()
!
      if (.not.allocated(sndbuf)) then
        allocate(sndbuf(mxmaglon+2,mxneed,nf),stat=ier) ! +2 for single halos
        if (ier /= 0)
     |    call shutdown('mp_mag_jslot: error allocating sndbuf')
      endif
      if (.not.allocated(rcvbuf)) then
        allocate(rcvbuf(mxmaglon+2,mxneed,nf),stat=ier)
        if (ier /= 0)
     |    call shutdown('mp_mag_jslot: error allocating rcvbuf')
      endif
      sndbuf = 0.
      rcvbuf = 0.
      njneed = 0
      do j=1,mxneed
        if (jneed(j) /= -1) njneed=njneed+1
      enddo
      if (any(jneed(1:njneed)==-1)) call shutdown('mp_mag_jslot jneed')
!
      call MPI_Comm_rank(cols_comm,tij,ier)
!     if (njneed >  0) then
!       write(6,"('mp_mag_jslot: tij=',i3,' nmagtaskj=',i3,' njneed=',
!    |    i3,' jneed=',/,(20i3))") tij,nmagtaskj,njneed,jneed(1:njneed)
!     endif ! njneed > 0
!
! Send needed lat indices to all tasks in my column:
! (redundant for alltoall)
      do n=0,nmagtaskj-1
        jneedall(:,n) = jneed(:)
      enddo

      call mpi_alltoall(jneedall,mxneed,MPI_INTEGER,
     |  peersneed,mxneed,MPI_INTEGER,cols_comm,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_mag_jslot call mpi_alltoall')
!
! Check if I have any needed lats, and who to send to:
      do n=0,nmagtaskj-1
        if (n==tij) cycle
        njhave = 0
        do j=1,mxneed
          if (peersneed(j,n) >= mlat00.and.peersneed(j,n) <= mlat11)then
            njhave = njhave+1 
            jhave(njhave) = peersneed(j,n)
            idest = n
            wid = itask_table(mytidi,idest)
          endif
        enddo
        if (njhave > 0) then
!         write(6,"('mp_mag_jslot: tij=',i3,' send to n=',i3,' njhave=',
!    |      i3,' jhave=',/,(20i4))") tij,n,njhave,jhave(1:njhave)
!
! Load send buffer:
          nlons = mlon11-mlon00+1
          do ifld=1,nf
            do j=1,njhave
              do i=mlon00,mlon11
                sndbuf(i-mlon00+1,j,ifld) = fin(i,jhave(j),ifld)
              enddo
!             write(6,"('mp_mag_jslot: tij=',i3,' Send jhave ',i4,
!    |          ' to task ',i3,' (world_comm task',i3,') ifld=',i3,
!    |          ' data=',/,(6e12.4))") tij,jhave(j),idest,wid,ifld,
!    |          sndbuf(1:nlons,j,ifld)
            enddo
          enddo
          len = nlons*njhave*nf
          msgid = mytid+wid

          call mpi_send(sndbuf(1:nlons,1:njhave,:),len,MPI_REAL8,
     |      idest,1,cols_comm,ier)

!         call mpi_isend(sndbuf(1:nlons,1:njhave,:),len,MPI_REAL8,
!    |      idest,1,cols_comm,msgid,ier)
!         call mpi_wait(mytid+idest,irstat,ier)
        endif
      enddo ! n=0,nmagtaskj-1
!
! Determine which tasks to receive which lats from. Task to
! receive from must be in same task column magtidi as I am.
      if (njneed > 0) then
        njhave = 0
        jhave(:) = -1
        do n=0,ntask-1
          njhave = 0
          do j=1,njneed
            if (jneed(j) >= tasks(n)%mlat0-1 .and. 
     |          jneed(j) <= tasks(n)%mlat1+1) then
              njhave = njhave+1
              jhave(njhave) = jneed(j)
            endif
          enddo
          if (njhave > 0 .and. tasks(n)%magtidi==magtidi) then
            isrc = tasks(n)%magtidj ! task id in cols_comm to recv from
            nlons = mlon11-mlon00+1
            len = nlons*njhave*nf
            msgid = mytid+n

            call mpi_recv(rcvbuf(1:nlons,1:njhave,:),len,MPI_REAL8,
     |        isrc,1,cols_comm,irstat,ier)

!           call mpi_irecv(rcvbuf(1:nlons,1:njhave,:),len,MPI_REAL8,
!    |        isrc,1,cols_comm,msgid,ier)
!           call mpi_wait(mytid+isrc,irstat,ier)
!
! Get data from receive buffer:
! real,intent(out) :: fout(mlon00:mlon11,mxneed) ! returned data at needed lats
            do ifld=1,nf
              do j=1,njhave
                nj = ixfind(jneed,mxneed,jhave(j),icount)
                if (nj==0) call shutdown('jhave(j) not in jneed') 
                do i=mlon00,mlon11
                  fout(i,nj,ifld) = rcvbuf(i-mlon00+1,j,ifld)
                enddo
!               write(6,"('mp_mag_jslot: tij=',i3,' Receive j=',i4,
!    |            ' jhave(j)=',i4,' from task ',i3,' (comm_world task ',
!    |            i3,') nj=',i4,' ifld=',i3,' data=',/,(6e12.4))") 
!    |            tij,j,jhave(j),isrc,n,nj,ifld,fout(:,nj,ifld)
              enddo ! j=1,njhave
            enddo ! ifld=1,nf
          endif ! jhave > 0
        enddo ! n=0,ntask-1
      endif ! njneed > 0
!
! Wait for completions:
!     ireq = (/isend,irecv/)
!     istat = 0
!     call mpi_waitall(2,ireq,istat,ier)
!     if (ier /= 0) call handle_mpi_err(ier,'mp_mag_jslot waitall')

      if (mpi_timing) then
        endtime = mpi_wtime()
        time_mag_jslot=time_mag_jslot+(endtime-starttime)
      endif
      end subroutine mp_mag_jslot
!-----------------------------------------------------------------------
      subroutine handle_mpi_err(ierrcode,string)
!
! Args:
      integer,intent(in) :: ierrcode
      character(len=*) :: string
!
! Local:
      character(len=80) :: errstring
      integer :: len_errstring,ier
!
      call mpi_error_string(ierrcode,errstring,len_errstring,ier)
      write(6,"(/,'>>> mpi error: ',a)") trim(string)
      write(6,"('  ierrcode=',i3,': ',a)") trim(errstring)
      end subroutine handle_mpi_err
!-----------------------------------------------------------------------
      subroutine report_mpi_timing
      use hist_module,only: nstep
      real :: time_totalmpi
!
! Called from main tgcm at end of run.
!
      if (.not.mpi_timing.or.nstep==0) return

      time_totalmpi =
     |  time_gather2root_prim+
     |  time_gather2root_sech+
     |  time_gather2root_f3d+
     |  time_gather2root_lbc+
     |  time_bndlats+
     |  time_bndlats_f2d+
     |  time_bndlons+
     |  time_bndlons_f3d+
     |  time_bndlons_f2d+
     |  time_polelats+
     |  time_polelat_f3d+
     |  time_gatherlons_f3d+
     |  time_scatterlons_f3d+
     |  time_periodic_f4d+
     |  time_periodic_f3d+
     |  time_periodic_f2d+
     |  time_bndlats_kmh+
     |  time_bndlons_kmh+
     |  time_mageq+
     |  time_mageq_jpm1+
     |  time_mageq_jpm3+
     |  time_magpole_2d+
     |  time_magpole_3d+
     |  time_magpoles+
     |  time_conjugate_points+
     |  time_mag_foldhem+
     |  time_mag_periodic_f2d+
     |  time_gather_pdyn+
     |  time_mag_halos+
     |  time_geo_halos+
     |  time_geo_halos_f3d+
     |  time_scatter_coeffs+
     |  time_scatter_phim

      write(6,"(/,72('-'))")

      write(6,"('Total wallclock time: mins=',f10.2,' hours=',f10.2)")
     |  time_totalrun/60.,time_totalrun/3600.

      write(6,"('Total MPI timing: mins=',e12.4,' hours=',e12.4,
     |  ' %Total runtime=',f9.2)")
     |  time_totalmpi/60.,time_totalmpi/3600.,
     |  time_totalmpi/time_totalrun*100.

      write(6,"('Subroutine         ',' Time (mins)',
     |  '  %Total mpi','  %Total run')")

      write(6,"('mp_gather2root_prim',e12.4,2f12.2)")
     | time_gather2root_prim/60.,
     | time_gather2root_prim/time_totalmpi*100.,
     | time_gather2root_prim/time_totalrun*100.

      write(6,"('mp_gather2root_sech',e12.4,2f12.2)")
     | time_gather2root_sech/60.,
     | time_gather2root_sech/time_totalmpi*100.,
     | time_gather2root_sech/time_totalrun*100.

      write(6,"('mp_gather2root_lbc ',e12.4,2f12.2)")
     | time_gather2root_lbc/60.,
     | time_gather2root_lbc/time_totalmpi*100.,
     | time_gather2root_lbc/time_totalrun*100.

      write(6,"('mp_bndlats         ',e12.4,2f12.2)")
     | time_bndlats/60.,
     | time_bndlats/time_totalmpi*100.,
     | time_bndlats/time_totalrun*100.

      write(6,"('mp_bndlats_f2d     ',e12.4,2f12.2)")
     | time_bndlats_f2d/60.,
     | time_bndlats_f2d/time_totalmpi*100.,
     | time_bndlats_f2d/time_totalrun*100.

      write(6,"('mp_bndlons         ',e12.4,2f12.2)")
     | time_bndlons/60.,
     | time_bndlons/time_totalmpi*100.,
     | time_bndlons/time_totalrun*100.

      write(6,"('mp_bndlons_f3d     ',e12.4,2f12.2)")
     | time_bndlons_f3d/60.,
     | time_bndlons_f3d/time_totalmpi*100.,
     | time_bndlons_f3d/time_totalrun*100.

      write(6,"('mp_polelats        ',e12.4,2f12.2)")
     | time_polelats/60.,
     | time_polelats/time_totalmpi*100.,
     | time_polelats/time_totalrun*100.

      write(6,"('mp_polelat_f3d     ',e12.4,2f12.2)")
     | time_polelat_f3d/60.,
     | time_polelat_f3d/time_totalmpi*100.,
     | time_polelat_f3d/time_totalrun*100.

      write(6,"('mp_gatherlons_f3d  ',e12.4,2f12.2)")
     | time_gatherlons_f3d/60.,
     | time_gatherlons_f3d/time_totalmpi*100.,
     | time_gatherlons_f3d/time_totalrun*100.

      write(6,"('mp_scatterlons_f3d ',e12.4,2f12.2)")
     | time_scatterlons_f3d/60.,
     | time_scatterlons_f3d/time_totalmpi*100.,
     | time_scatterlons_f3d/time_totalrun*100.

      write(6,"('mp_periodic_f4d    ',e12.4,2f12.2)")
     | time_periodic_f4d/60.,
     | time_periodic_f4d/time_totalmpi*100.,
     | time_periodic_f4d/time_totalrun*100.

      write(6,"('mp_periodic_f3d    ',e12.4,2f12.2)")
     | time_periodic_f3d/60.,
     | time_periodic_f3d/time_totalmpi*100.,
     | time_periodic_f3d/time_totalrun*100.

      write(6,"('mp_periodic_f2d    ',e12.4,2f12.2)")
     | time_periodic_f2d/60.,
     | time_periodic_f2d/time_totalmpi*100.,
     | time_periodic_f2d/time_totalrun*100.

      write(6,"('mp_bndlats_kmh     ',e12.4,2f12.2)")
     | time_bndlats_kmh/60.,
     | time_bndlats_kmh/time_totalmpi*100.,
     | time_bndlats_kmh/time_totalrun*100.

      write(6,"('mp_bndlons_kmh     ',e12.4,2f12.2)")
     | time_bndlons_kmh/60.,
     | time_bndlons_kmh/time_totalmpi*100.,
     | time_bndlons_kmh/time_totalrun*100.

      write(6,"('mp_mageq           ',e12.4,2f12.2)")
     | time_mageq/60.,
     | time_mageq/time_totalmpi*100.,
     | time_mageq/time_totalrun*100.
      write(6,"('mp_mageq_jpm1      ',e12.4,2f12.2)")
     | time_mageq_jpm1/60.,
     | time_mageq_jpm1/time_totalmpi*100.,
     | time_mageq_jpm1/time_totalrun*100.

      write(6,"('mp_mageq_jpm3      ',e12.4,2f12.2)")
     | time_mageq_jpm3/60.,
     | time_mageq_jpm3/time_totalmpi*100.,
     | time_mageq_jpm3/time_totalrun*100.

      write(6,"('mp_magpole_2d      ',e12.4,2f12.2)")
     | time_magpole_2d/60.,
     | time_magpole_2d/time_totalmpi*100.,
     | time_magpole_2d/time_totalrun*100.

      write(6,"('mp_magpole_3d      ',e12.4,2f12.2)")
     | time_magpole_3d/60.,
     | time_magpole_3d/time_totalmpi*100.,
     | time_magpole_3d/time_totalrun*100.

      write(6,"('mp_magpoles        ',e12.4,2f12.2)")
     | time_magpoles/60.,
     | time_magpoles/time_totalmpi*100.,
     | time_magpoles/time_totalrun*100.

      write(6,"('mp_conjugate_points',e12.4,2f12.2)")
     | time_conjugate_points/60.,
     | time_conjugate_points/time_totalmpi*100.,
     | time_conjugate_points/time_totalrun*100.

      write(6,"('mp_foldhem         ',e12.4,2f12.2)")
     | time_mag_foldhem/60.,
     | time_mag_foldhem/time_totalmpi*100.,
     | time_mag_foldhem/time_totalrun*100.

      write(6,"('mp_mag_periodic_f2d',e12.4,2f12.2)")
     | time_mag_periodic_f2d/60.,
     | time_mag_periodic_f2d/time_totalmpi*100.,
     | time_mag_periodic_f2d/time_totalrun*100.

      write(6,"('mp_gather_pdyn     ',e12.4,2f12.2)")
     | time_gather_pdyn/60.,
     | time_gather_pdyn/time_totalmpi*100.,
     | time_gather_pdyn/time_totalrun*100.

      write(6,"('mp_mag_halos       ',e12.4,2f12.2)")
     | time_mag_halos/60.,
     | time_mag_halos/time_totalmpi*100.,
     | time_mag_halos/time_totalrun*100.

      write(6,"('mp_geo_halos       ',e12.4,2f12.2)")
     | time_geo_halos/60.,
     | time_geo_halos/time_totalmpi*100.,
     | time_geo_halos/time_totalrun*100.

      write(6,"('mp_geo_halos_f3d   ',e12.4,2f12.2)")
     | time_geo_halos_f3d/60.,
     | time_geo_halos_f3d/time_totalmpi*100.,
     | time_geo_halos_f3d/time_totalrun*100.

      write(6,"('mp_scatter_coeffs  ',e12.4,2f12.2)")
     | time_scatter_coeffs/60.,
     | time_scatter_coeffs/time_totalmpi*100.,
     | time_scatter_coeffs/time_totalrun*100.

      write(6,"('mp_scatter_phim    ',e12.4,2f12.2)")
     | time_scatter_phim/60.,
     | time_scatter_phim/time_totalmpi*100.,
     | time_scatter_phim/time_totalrun*100.

      write(6,"('mp_gather_f2d     ',e12.4,2f12.2)")
     | time_gather_f2d/60.,
     | time_gather_f2d/time_totalmpi*100.,
     | time_gather_f2d/time_totalrun*100.

      write(6,"('mp_scatter_f2d     ',e12.4,2f12.2)")
     | time_scatter_f2d/60.,
     | time_scatter_f2d/time_totalmpi*100.,
     | time_scatter_f2d/time_totalrun*100.

      write(6,"(/,72('-'))")
      end subroutine report_mpi_timing
!-----------------------------------------------------------------------
#endif
      end module mpi_module
