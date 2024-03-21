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
      use params_module,only: nlevp1
      use fields_module,only: nf4d,foutput,fzg,fnrh,fkmh
      use mpi_module,only: lon0,lon1,lat0,lat1
      use chemrates_module,only: alloc_tdep
      use qrj_module,only: alloc_q
      use n2d_module,only: alloc_n2d
      use n4s_module,only: alloc_n4s
      use no_module,only: alloc_no
      use ar_module,only: alloc_ar
      use magpres_g,only: alloc_pg
      use gswm_module,only: alloc_gswm
      use bgrd_data_module,only: alloc_bgrddata
      use ctmt_module,only: alloc_ctmt
      use pdynamo_module,only: alloc_pdyn
      use dyndiag_module,only: alloc_dyndiag
      implicit none
!
! Local:
      integer :: istat
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
      call alloc_n2d(lon0,lon1,lat0,lat1)
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
      call alloc_dyndiag(lon0-2,lon1+2,lat0-2,lat1+2)
!
! We are doing parallel i/o to netcdf output histories,
! so when it comes time to write a history,
! the task subdomains are directly written to the file.
!
! The pointer is declared in fields.F:
!     real,pointer :: foutput(:,:,:,:) ! (nlevp1,lon0:lon1,lat0:lat1,nf4d)
! 
      allocate(foutput(nlevp1,lon0:lon1,lat0:lat1,nf4d),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating 3d data for ',
     |    'all fields on the current task: istat=',i3)") istat
        call shutdown('allocdata')
      endif
!
      allocate(fzg(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating fzg')")
        call shutdown('allocdata')
      endif
!
      allocate(fnrh(nlevp1,lon0-2:lon1+2,lat0-2:lat1+2),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating fnrh')")
        call shutdown('allocdata')
      endif
!
      allocate(fkmh(nlevp1,lon0-2:lon1+2,lat0-2:lat1+2),stat=istat)
      if (istat /= 0) then
        write(6,"('>>> allocdata: error allocating fkmh')")
        call shutdown('allocdata')
      endif
!
      end subroutine allocdata
