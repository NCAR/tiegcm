!
      subroutine allocdata
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Allocate msc data. Main fields are not allocated here, see
! init_fields in fields.F. This is not in field.F in order to 
! avoid circular module dependencies.
!
      use params_module,only: nlevp1,nlonp4,nlat
      use fields_module,only: f4d,f3d,nf4d,nf3d,field_initval,foutput,
     |  fzg
      use mpi_module,only: lon0,lon1,lat0,lat1,mytid
      use chemrates_module,only: alloc_tdep
      use qrj_module,only: alloc_q
      use n4s_module,only: alloc_n4s
      use no_module,only: alloc_no
      use ar_module,only: alloc_ar
      use magpres_g,only: alloc_pg
      use gswm_module,only: alloc_gswm
      use bgrd_data_module,only: alloc_bgrddata
      use ctmt_module,only: alloc_ctmt
      use pdynamo_module,only: alloc_pdyn
      implicit none
!
! Local:
      integer :: i,istat,nlons,nlats,len 
!
! Allocate temperature-dependent chemical reaction rates:
      call alloc_tdep
!
! Allocate ionization rates (qrj module):
      call alloc_q(lon0,lon1,lat0,lat1)
!
! Allocate J(mag.pressure,gravity (magpres_g module):
      call alloc_pg(lon0,lon1,lat0,lat1)
!
! Allocate boundary conditions and production and loss terms
! for minor composition species:
      call alloc_n4s(lon0,lon1,lat0,lat1)
      call alloc_no (lon0,lon1,lat0,lat1)
      call alloc_ar (lon0,lon1,lat0,lat1)
!
! Allocate for gswm:
      call alloc_gswm(lon0,lon1,lat0,lat1)
!
! Allocate for bgrddata:
      call alloc_bgrddata(lon0,lon1,lat0,lat1)
!
! Allocate for ctmt:
      call alloc_ctmt(lon0,lon1,lat0,lat1)
!
! Allocate needed fields for parallel dynamo (pdynamo.F):
      call alloc_pdyn
!
! We are not doing parallel i/o to netcdf output histories, so when
! it comes time to write a history, the task subdomains are gathered
! to the root task (see call mp_gather2root in advance), and only
! the root task writes to the file. Here we allocate the full 3d grid
! for all history fields on the root task. This is where the data
! will be gathered and written by the root when a history is written. 
!
! The pointer is declared in fields.F:
!     real,pointer :: foutput(:,:,:,:) ! (nlevp1,nlonp4,nlat,nf4d)
! 
#ifdef MPI
      if (mytid==0) then
        allocate(foutput(nlevp1,nlonp4,nlat,nf4d),stat=istat)
        if (istat /= 0) then
          write(6,"('>>> allocdata: error allocating 3d data for ',
     |      'all fields on root task: nlevp1=',i3,' nlonp4=',i3,
     |      ' nlat=',i3,' nf4d=',i3,' istat=',i3)")
     |      nlevp1,nlonp4,nlat,nf4d,istat
        else
          write(6,"('Allocated 3d data for all fields on root task:')")
          len = nlevp1*nlonp4*nlat*nf4d
          write(6,"('  nlevp1*nlonp4*nlat*nf4d = ',i3,'*',i3,'*',
     |      i3,'*',i3,' = ',i10,' words *8 = ',i12,' bytes.')")
     |      nlevp1,nlonp4,nlat,nf4d,len,len*8
        endif
      endif
      allocate(fzg(nlevp1,nlonp4,nlat),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating fzg')")
      else
        write(6,"('allocdata: all tasks allocate fzg: nlevp1=',i4,
     |    ' nlonp4=',i4,' nlat=',i4)") nlevp1,nlonp4,nlat
      endif
#else
      allocate(foutput(nlevp1,nlonp4,nlat,nf4d),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating 3d data for ',
     |    'all fields on root task: nlevp1=',i3,' nlonp4=',i3,
     |    ' nlat=',i3,' nf4d=',i3,' istat=',i3)")
     |    nlevp1,nlonp4,nlat,nf4d,istat
      else
        write(6,"('Allocated 3d data for all fields on root task:')")
        len = nlevp1*nlonp4*nlat*nf4d
        write(6,"('  nlevp1*nlonp4*nlat*nf4d = ',i3,'*',i3,'*',
     |    i3,'*',i3,' = ',i10,' words *8 = ',i12,' bytes.')")
     |    nlevp1,nlonp4,nlat,nf4d,len,len*8
      endif
      allocate(fzg(nlevp1,nlonp4,nlat),stat=istat)
      if (istat /= 0) 
     |  write(6,"('>>> allocdata: error allocating fzg')")
#endif
      end subroutine allocdata
