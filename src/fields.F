!
      module fields_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! There are allocatable arrays and an array of field structures for 
!   3-d and 4-d fields. Subdomains are allocated to the allocatable 
!   arrays, and structure data pointers point to the allocatable arrays,
!   e.g. f4d(i_tn)%data(:,:,:,:) => tn(:,:,:,:), where the 4 dimensions
!   are (pressure,longitude,latitude,2). The last dimension is for
!   previous and current time steps. 3-d fields are the same as 4-d
!   without the final time dimension. 
!
      use params_module,only: nlonp1,nlonp4,nlat,nlatp1,nlevp1,
     |  mxfsech,spval,nmlat,nmlonp1,mlev0,mlev1
      implicit none
      integer,parameter :: 
     |  nf4d          = 30, ! number of 4-d fields
     |  nf4d_hist     = 28, ! number of 4-d fields on primary histories
     |  nf3d          = 46  ! number of 3-d fields
      real,parameter ::
     |  field_initval = 0.  ! initialization value for fields data
      integer,parameter :: 
     |  longname_len  = 80, ! length of field long name 
     |  shortname_len = 16, ! length of field short name 
     |  units_len     = 16  ! length of field units attribute
      integer :: levd0,levd1, lond0,lond1, latd0,latd1 ! subdomain dimensions
!
! 4-d fields structure type:
! Data pointers will point to allocatable target arrays.
!
      type fields_4d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        character(len=16)            :: vcoord ! midpoints or interfaces
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi          ! flag used by some mpi routines
        real                         :: polesign
        real,pointer                 :: data(:,:,:,:)
      end type fields_4d
      type (fields_4d) :: f4d(nf4d)
!
! Indices, to f4d fields, e.g. f4d(i_tn)%data is neutral temperature.
      integer ::
     |  i_tn    ,i_un     ,i_vn     ,i_o2     ,i_o1     ,
     |  i_n4s   ,i_no     ,i_op     ,i_n2d    ,i_ti     ,
     |  i_te    ,i_ne     ,i_o2p    ,i_w      ,i_z      ,
     |  i_poten ,i_tn_nm  ,i_un_nm  ,i_vn_nm  ,i_o2_nm  ,
     |  i_o1_nm ,i_n4s_nm ,i_no_nm  ,i_op_nm  ,i_barm   ,
     |  i_vc    ,i_he     ,i_he_nm  ,i_ar     ,i_ar_nm
!
! Allocatable target arrays for 4-d fields:
! Allocation will be:   tn(levd0:levd1,lond0 :lond1 ,latd0 :latd1 ,2)
! which is the same as: tn(1    :nlev ,lon0-2:lon1+2,lat0-2:lat1+2,2)
!
      real,dimension(:,:,:,:),target,allocatable :: ! (k,i,j,2)
     |  tn      ,un       ,vn       ,o2       ,o1       ,
     |  n4s     ,no       ,op       ,n2d      ,ti       ,
     |  te      ,ne       ,o2p      ,w        ,z        ,
     |  poten   ,tn_nm    ,un_nm    ,vn_nm    ,o2_nm    ,
     |  o1_nm   ,n4s_nm   ,no_nm    ,op_nm    ,barm     ,
     |  vc      ,he       ,he_nm    ,ar       ,ar_nm
!
! 3-d fields structure type:
! Data pointers will point to allocatable target arrays.
!
      type fields_3d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        character(len=16)            :: vcoord ! midpoints or interfaces
        character(len=8)             :: dimnames(3) ! e.g., lat,lon,lev, or lat,lon,' '
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi,         ! flag used by some mpi routines
     |    task0_only   ! if true, field was defined at root task only
        real                         :: polesign
        real,pointer                 :: data(:,:,:)
      end type fields_3d
      type (fields_3d) :: f3d(nf3d)
!
! 3-d fields and indices, e.g., reference f3d(kldt)%data
      integer :: 
     |  i_kldt   ,i_kldu    ,i_kldv    ,i_kldo2 ,i_kldo1 ,i_kldhe,
     |  i_xnmbar ,i_xnmbari ,i_xnmbarm ,i_cp    ,i_kt    ,
     |  i_km     ,i_ui      ,i_vi      ,i_wi    ,i_vo2   ,
     |  i_vo1    ,i_vn2     ,i_sco2    ,i_sco1  ,i_scn2  ,
     |  i_xiop2p ,i_xiop2d  ,i_nplus   ,i_n2p   ,i_nop   ,
     |  i_lxx    ,i_lyy     ,i_lxy     ,i_lyx   ,i_qji_ti,
     |  i_qji_tn ,i_cool_implicit      ,i_cool_explicit  ,
     |  i_hdt    ,i_hdu     ,i_hdv     ,i_hdo2  ,i_hdo1  ,i_hdhe,
     |  i_ped    ,i_hall    ,i_lam1    ,i_zg    ,i_n2    ,i_wn
!
! Allocatable target arrays 3-d:
! Allocation will be:   tn(levd0:levd1,lond0 :lond1 ,latd0 :latd1 )
! which is the same as: tn(1    :nlev ,lon0-2:lon1+2,lat0-2:lat1+2)
!
      real,dimension(:,:,:),target,allocatable :: ! (k,i,j)
     |  kldt     ,kldu     ,kldv     ,kldo2    ,kldo1   ,kldhe,
     |  xnmbar   ,xnmbari  ,xnmbarm  ,cp       ,kt      ,
     |  km       ,ui       ,vi       ,wi       ,vo2     ,
     |  vo1      ,vn2      ,sco2     ,sco1     ,scn2    ,
     |  xiop2p   ,xiop2d   ,nplus    ,n2p      ,nop     ,
     |  lxx      ,lyy      ,lxy      ,lyx      ,qji_ti  ,
     |  qji_tn   ,cool_implicit      ,cool_explicit     ,
     |  hdt      ,hdu      ,hdv      ,hdo2     ,hdo1    ,hdhe,
     |  ped      ,hall     ,lam1     ,zg       ,n2      ,wn

!
! 2-d field type (used for 2d secondary history fields):
      type fields_2d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    prognostic,  ! true if prognostic (diagnostic if false)
     |    mpi,         ! flag used by some mpi routines
     |    task0_only   ! if true, field was defined at root task only
        real                         :: polesign
        real,pointer                 :: data(:,:) ! (k,i)
      end type fields_2d
!
! Electric potential on geographic and magnetic grids:
! (full domains until dynamo is parallelized)
      real ::
     |  phim3d(nmlonp1,nmlat,mlev0:mlev1),  ! 3d electric potential magnetic
     |  emphi3d(nmlonp1,nmlat,mlev0:mlev1), ! 3d eastward electric field magnetic
     |  emlam3d(nmlonp1,nmlat,mlev0:mlev1), ! 3d equatorw. electric field magnetic
     |  emz3d(nmlonp1,nmlat,mlev0:mlev1)    ! 3d upward (?) electric field magnetic
!
! Secondary history field structures.
! These are initialized in sub init_fsech, and data is set in addfld.F.
      type fields_sech
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: task0_only
        logical                      :: prognostic
        logical                      :: mag,geo
        integer                      :: ndims        ! 2d or 3d
        character(len=8)             :: dimnames(3)  ! dims on history
        integer                      :: dimsizes(3)  ! dim sizes on history
        real,pointer                 :: data(:,:,:)  ! allocated by addfld
      end type fields_sech
      type(fields_sech) :: fsechist(mxfsech)
!
! Time indices for rightmost dimension of 4d data, itp for current timestep, 
! itc for next timestep. Fields at the previous timestep (time n-1) are
! saved at both itp and itc (e.g., tn_nm, un_nm, etc).
!
      integer :: itc,itp 
!
! If fakeflds is true, use fake dimensions for fields
! (for dry runs, testing, etc)
!
      logical,parameter :: fakeflds = .false.
!
      real,dimension(nlevp1,nlonp4+1,-2:nlat) ::
     |  fnrh,  ! eddy viscosity
     |  fkmh   ! M/T
!
! Lower boundary conditions (bottom interface level) for t,u,v from 
! source history, and at current model time:
      real,allocatable,dimension(:,:) ::
     |  tlbc,     ulbc,     vlbc,     ! subdomains (lond0:lond1,latd0:latd1)
     |  tlbc_nm,  ulbc_nm,  vlbc_nm,  ! subdomains (lond0:lond1,latd0:latd1)
     |  tlbc_glb, ulbc_glb, vlbc_glb, ! global (nlonp4,nlat)
     |  tlbc_nm_glb, ulbc_nm_glb, vlbc_nm_glb  ! global (nlonp4,nlat)
!
! Full 3d grid with all primary history fields for writing to netcdf
! history files. This will be allocated only on the root task 
! (see allocdata):
      real,allocatable :: foutput(:,:,:,:) ! (nlevp1,nlonp4,nlat,nf4d_hist)
      real,allocatable :: fzg(:,:,:)       ! (nlevp1,nlonp4,nlat)
!
! Magnetic field diagnostics (see getapex.F, apex.F90, diags.F)
!     real,dimension(nlonp4,nlat) ::
!    |  bx,      ! Magnetic field eastward component (Gauss)
!    |  by,      ! Magnetic field northward component (Gauss)
!    |  bz,      ! Magnetic field upward component (Gauss)
!    |  bmag     ! Magnitude of magnetic field (Gauss)

      contains
!-----------------------------------------------------------------------
!
      subroutine init_4d(lon0,lon1,lat0,lat1,mytid,iprint)
!
! Set names, units, indices and pointers for f4d and f3d field structures,
! and allocate 3d and 4d field arrays. Also make other data allocations.
!
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,mytid,iprint
!
! Local:
      integer :: n,istat,i
!
! Fields are allocated at full task subdomain, including ghost cells:
      levd0 = 1      ; levd1 = nlevp1
      lond0 = lon0-2 ; lond1 = lon1+2
      latd0 = lat0-2 ; latd1 = lat1+2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      n = 1
!
! 4-d fields (long and short names, units, index):
! (pointer definition must follow allocate statement)
! Note it is difficult to make a subroutine to do the allocations because
!   allocatable arrays cannot be dummy arguments. (Once they are allocated,
!   they can be passed to subroutines and declared as real, but subroutine 
!   dummy arguments cannot be declared as allocatable arrays.)
!
      allocate(tn(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL TEMPERATURE')")
      f4d(n)%short_name = "TN"
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => tn
      i_tn = n ; n = n+1
!
      allocate(un(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND (+EAST)')")
      f4d(n)%short_name = "UN"
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => un
      i_un = n ; n = n+1
!
      allocate(vn(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND (+NORTH)')")
      f4d(n)%short_name = "VN"
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vn
      i_vn = n ; n = n+1
!
      allocate(o2(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN')")
      f4d(n)%short_name = ("O2")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o2
      i_o2 = n ; n = n+1
!
      allocate(o1(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN')")
      f4d(n)%short_name = ("O1")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o1
      i_o1 = n ; n = n+1
!
      allocate(n4s(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S')")
      f4d(n)%short_name = ("N4S")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n4s
      i_n4s = n ; n = n+1
!
      allocate(no(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NITRIC OXIDE')")
      f4d(n)%short_name = ("NO")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => no
      i_no = n ; n = n+1
!
      allocate(he(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('HELIUM')")
      f4d(n)%short_name = ("HE")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => he
      i_he = n ; n = n+1
!
      allocate(ar(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ARGON (AR)')")
      f4d(n)%short_name = ("AR")
      f4d(n)%units = "MMR"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => ar
      i_ar = n ; n = n+1
!
      allocate(op(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O+ ION')")
      f4d(n)%short_name = ("OP")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => op
      i_op = n ; n = n+1
!
      allocate(n2d(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N2D')")
      f4d(n)%short_name = ("N2D")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n2d
      i_n2d = n ; n = n+1
!
      allocate(ti(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ION TEMPERATURE')")
      f4d(n)%short_name = ("TI")
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => ti
      i_ti = n ; n = n+1
!
      allocate(te(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRON TEMPERATURE')")
      f4d(n)%short_name = ("TE")
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => te
      i_te = n ; n = n+1
!
      allocate(ne(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRON DENSITY')")
      f4d(n)%short_name = ("NE")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => ne
      i_ne = n ; n = n+1
!
      allocate(o2p(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O2+ ION')")
      f4d(n)%short_name = ("O2P")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o2p
      i_o2p = n ; n = n+1
!
! 10/27/05 btf: "W" is replaced by "OMEGA":
!
!     allocate(w(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
!     write(f4d(n)%long_name,"('VERTICAL VELOCITY (PLUS UP)')")
!     f4d(n)%short_name = ("W")
!     f4d(n)%units = "cm/s"
!     f4d(n)%vcoord = "interfaces"
!     f4d(n)%data => w
!     i_w = n ; n = n+1
!
      allocate(w(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('VERTICAL MOTION')")
      f4d(n)%short_name = ("OMEGA")
      f4d(n)%units = "s-1"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => w
      i_w = n ; n = n+1
!
      allocate(z(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('GEOPOTENTIAL HEIGHT')")
      f4d(n)%short_name = ("Z")
      f4d(n)%units = "cm"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => z
      i_z = n ; n = n+1
!
      allocate(poten(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRIC POTENTIAL')")
      f4d(n)%short_name = ("POTEN")
      f4d(n)%units = "volts"
      f4d(n)%vcoord = "interfaces"
      f4d(n)%data => poten
      i_poten = n ; n = n+1
!
      allocate(tn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL TEMPERATURE (TIME N-1)')")
      f4d(n)%short_name = ("TN_NM")
      f4d(n)%units = "K"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => tn_nm
      i_tn_nm = n ; n = n+1
!
      allocate(un_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("UN_NM")
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => un_nm
      i_un_nm = n ; n = n+1
!
      allocate(vn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("VN_NM")
      f4d(n)%units = "cm/s"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vn_nm
      i_vn_nm = n ; n = n+1
!
      allocate(o2_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O2_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o2_nm
      i_o2_nm = n ; n = n+1
!
      allocate(o1_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O1_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => o1_nm
      i_o1_nm = n ; n = n+1
!
      allocate(n4s_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S (TIME N-1)')")
      f4d(n)%short_name = ("N4S_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => n4s_nm
      i_n4s_nm = n ; n = n+1
!
      allocate(no_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NO (TIME N-1)')")
      f4d(n)%short_name = ("NO_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => no_nm
      i_no_nm = n ; n = n+1
!
      allocate(op_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('OP (TIME N-1)')")
      f4d(n)%short_name = ("OP_NM")
      f4d(n)%units = "cm-3"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => op_nm
      i_op_nm = n ; n = n+1
!
      allocate(he_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('HELIUM (TIME N-1)')")
      f4d(n)%short_name = ("HE_NM")
      f4d(n)%units = "mmr"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => he_nm
      i_he_nm = n ; n = n+1
!
      allocate(ar_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ARGON (TIME N-1)')")
      f4d(n)%short_name = ("AR_NM")
      f4d(n)%units = "MMR"
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => ar_nm
      i_ar_nm = n ; n = n+1
!
      allocate(barm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MEAN MOLECULAR WEIGHT')")
      f4d(n)%short_name = ("BARM")
      f4d(n)%units = ' '
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => barm
      i_barm = n ; n = n+1
!
      allocate(vc(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('COS(PHI)*VN')")
      f4d(n)%short_name = ("VC")
      f4d(n)%units = ' '
      f4d(n)%vcoord = "midpoints"
      f4d(n)%data => vc
      i_vc = n ; n = n+1
!
! Check number of fields allocated:
      if (n-1 /= nf4d) then
        write(6,"(/,'>>> init_4d: wrong number of 4-d fields?',
     |    ' n-1=',i5,' nf4d=',i5)") n-1,nf4d
      endif
!
! itp,itc are indices to rightmost dimension of field data, for
! previous and current time steps.
      itp = 1
      itc = 2
!
! Set polesign for crossing the poles (see mk_polelats and mp_bndlats)
!
! 11/05 btf: Getting the following error here, attempting to initialize
!            f4d(:)%polesign, from hao (callisto) pgf90 5.1-6:
!
! Lowering Error: symbol data$sd is a member reference
! Lowering Error: symbol data$sd is a member reference
! PGF90-F-0000-Internal compiler error. Errors in Lowering      21 (/home/foster/tgcm/tiegcm-hist-lbc/src/fields.F: 437)
! PGF90/any Linux/x86 5.1-6: compilation aborted
!
! However, this error does not occur if fields.F is local!
!
!     f4d(:)%polesign = 1. ! 11/05 btf: pgf90 at hao does not like this
!     do i=1,size(f4d)
      do i=1,nf4d
        f4d(i)%polesign = 1.
      enddo
      f4d(i_un)%polesign = -1.
      f4d(i_un_nm)%polesign = -1.
      f4d(i_vn)%polesign = -1.
      f4d(i_vn_nm)%polesign = -1.
      f4d(i_n2d)%polesign = 0.
      f4d(i_ne )%polesign = 0.
      f4d(i_o2p)%polesign = 0.
!
! All f4d fields are on geographic grid:
      f4d%magnetic = .false.
!
! f4d fields are "prognostic":
      f4d%prognostic = .true.
!
! Init data:
      do i=1,nf4d
        f4d(i)%data = field_initval
      enddo
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 4-d fields:')") nf4d
        do n=1,nf4d
          call print_f4d(f4d(n),n)
        enddo
      endif
      end subroutine init_4d
!-----------------------------------------------------------------------
      subroutine init_3d(lon0,lon1,lat0,lat1,mytid,iprint)
!
! Set names, units, indices and pointers for f4d and f3d field structures,
! and allocate 3d and 4d field arrays. Also make other data allocations.
!
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1,mytid,iprint
!
! Local:
      integer :: n,istat
!
      n = 1
!
! 3-d fields (long and short names, units, index):
!
      allocate(kldt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('TN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDT"
      f3d(n)%units = " "
      f3d(n)%data => kldt
      i_kldt = n ; n = n+1
!
      allocate(kldu(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('UN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDU"
      f3d(n)%units = " "
      f3d(n)%data => kldu
      i_kldu = n ; n = n+1
!
      allocate(kldv(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('VN HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDV"
      f3d(n)%units = " "
      f3d(n)%data => kldv
      i_kldv = n ; n = n+1
!
      allocate(kldo2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDO2"
      f3d(n)%units = " "
      f3d(n)%data => kldo2
      i_kldo2 = n ; n = n+1
!
      allocate(kldo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDO1"
      f3d(n)%units = " "
      f3d(n)%data => kldo1
      i_kldo1 = n ; n = n+1
!
      allocate(kldhe(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('HE HORIZONTAL DIFFUSION COEFF')")
      f3d(n)%short_name = "KLDHE"
      f3d(n)%units = " "
      f3d(n)%data => kldhe
      i_kldhe = n ; n = n+1

      allocate(xnmbar(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('p0*e(-z)/kT*barm')")
      f3d(n)%short_name = ("XNMBAR")
      f3d(n)%units = ' '
      f3d(n)%data => xnmbar
      i_xnmbar = n ; n = n+1
!
      allocate(xnmbari(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('p0*e(-z)/kT*barm at interfaces')")
      f3d(n)%short_name = ("XNMBARI")
      f3d(n)%units = ' '
      f3d(n)%data => xnmbari
      i_xnmbari = n ; n = n+1
!
      allocate(xnmbarm(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('p0*e(-z)/kT*barm at midpoints')")
      f3d(n)%short_name = ("XNMBARM")
      f3d(n)%units = ' '
      f3d(n)%data => xnmbarm
      i_xnmbarm = n ; n = n+1
!
      allocate(cp(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('SPECIFIC HEAT')")
      f3d(n)%short_name = "CP"
      f3d(n)%units = " "
      f3d(n)%data => cp
      i_cp = n ; n = n+1
!
      allocate(kt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MOLECULAR THERMAL CONDUCTIVITY')")
      f3d(n)%short_name = "KT"
      f3d(n)%units = " "
      f3d(n)%data => kt
      i_kt = n ; n = n+1
!
      allocate(km(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MOLECULAR DIFFUSION')")
      f3d(n)%short_name = "KM"
      f3d(n)%units = " "
      f3d(n)%data => km
      i_km = n ; n = n+1
!
      allocate(ui(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('ZONAL ION DRIFT VELOCITY')")
      f3d(n)%short_name = "UI"
      f3d(n)%units = " "
      f3d(n)%data => ui
      i_ui = n ; n = n+1
!
      allocate(vi(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('MERIDIONAL ION DRIFT VELOCITY')")
      f3d(n)%short_name = "VI"
      f3d(n)%units = " "
      f3d(n)%data => vi
      i_vi = n ; n = n+1
!
      allocate(wi(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('VERTICAL ION DRIFT VELOCITY')")
      f3d(n)%short_name = "WI"
      f3d(n)%units = " "
      f3d(n)%data => wi
      i_wi = n ; n = n+1
!
      allocate(vo2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 LINE INTEGRAL')")
      f3d(n)%short_name = "VO2"
      f3d(n)%units = " "
      f3d(n)%data => vo2
      i_vo2 = n ; n = n+1
!
      allocate(vo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 LINE INTEGRAL')")
      f3d(n)%short_name = "VO1"
      f3d(n)%units = " "
      f3d(n)%data => vo1
      i_vo1 = n ; n = n+1
!
      allocate(vn2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2 LINE INTEGRAL')")
      f3d(n)%short_name = "VN2"
      f3d(n)%units = " "
      f3d(n)%data => vn2
      i_vn2 = n ; n = n+1
!
      allocate(sco2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 VERTICAL COLUMN DENSITY')")
      f3d(n)%short_name = "SCO2"
      f3d(n)%units = " "
      f3d(n)%data => sco2
      i_sco2 = n ; n = n+1
!
      allocate(sco1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 VERTICAL COLUMN DENSITY')")
      f3d(n)%short_name = "SCO1"
      f3d(n)%units = " "
      f3d(n)%data => sco1
      i_sco1 = n ; n = n+1
!
      allocate(scn2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2 VERTICAL COLUMN DENSITY')")
      f3d(n)%short_name = "SCN2"
      f3d(n)%units = " "
      f3d(n)%data => scn2
      i_scn2 = n ; n = n+1
!
      allocate(xiop2p(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('XIOP2P')")
      f3d(n)%short_name = "XIOP2P"
      f3d(n)%units = " "
      f3d(n)%data => xiop2p
      i_xiop2p = n ; n = n+1
!
      allocate(xiop2d(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('XIOP2D')")
      f3d(n)%short_name = "XIOP2D"
      f3d(n)%units = " "
      f3d(n)%data => xiop2d
      i_xiop2d = n ; n = n+1
!
      allocate(nplus(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N+ ION')")
      f3d(n)%short_name = "NPLUS"
      f3d(n)%units = " "
      f3d(n)%data => nplus
      i_nplus = n ; n = n+1
!
      allocate(n2p(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('N2+ ION')")
      f3d(n)%short_name = "N2P"
      f3d(n)%units = " "
      f3d(n)%data => n2p
      i_n2p = n ; n = n+1
!
      allocate(nop(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('NO+ ION')")
      f3d(n)%short_name = "NOP"
      f3d(n)%units = " "
      f3d(n)%data => nop
      i_nop = n ; n = n+1
!
      allocate(lxx(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('LAMDA ION DRAG XX')")
      f3d(n)%short_name = "LXX"
      f3d(n)%units = " "
      f3d(n)%data => lxx
      i_lxx = n ; n = n+1
!
      allocate(lyy(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('LAMDA ION DRAG YY')")
      f3d(n)%short_name = "LYY"
      f3d(n)%units = " "
      f3d(n)%data => lyy
      i_lyy = n ; n = n+1
!
      allocate(lxy(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('LAMDA ION DRAG XY')")
      f3d(n)%short_name = "LXY"
      f3d(n)%units = " "
      f3d(n)%data => lxy
      i_lxy = n ; n = n+1
!
      allocate(lyx(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('LAMDA ION DRAG YX')")
      f3d(n)%short_name = "LYX"
      f3d(n)%units = " "
      f3d(n)%data => lyx
      i_lyx = n ; n = n+1
!
      allocate(qji_ti(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('ION JOULE HEATING FOR TI')")
      f3d(n)%short_name = "QJI_TI"
      f3d(n)%units = " "
      f3d(n)%data => qji_ti
      i_qji_ti = n ; n = n+1
!
      allocate(qji_tn(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('ION JOULE HEATING FOR TN')")
      f3d(n)%short_name = "QJI_TN"
      f3d(n)%units = " "
      f3d(n)%data => qji_tn
      i_qji_tn = n ; n = n+1
!
      allocate(cool_implicit(levd0:levd1,lond0:lond1,latd0:latd1),
     |  stat=istat)
      write(f3d(n)%long_name,"('IMPLICIT COOLING TERM')")
      f3d(n)%short_name = "cool_implicit"
      f3d(n)%units = " "
      f3d(n)%data => cool_implicit
      i_cool_implicit = n ; n = n+1
!
      allocate(cool_explicit(levd0:levd1,lond0:lond1,latd0:latd1),
     |  stat=istat)
      write(f3d(n)%long_name,"('EXPLICIT COOLING TERM')")
      f3d(n)%short_name = "cool_explicit"
      f3d(n)%units = " "
      f3d(n)%data => cool_explicit
      i_cool_explicit = n ; n = n+1
!
      allocate(hdt(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('TN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdt"
      f3d(n)%units = " "
      f3d(n)%data => hdt
      i_hdt = n ; n = n+1
!
      allocate(hdu(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('UN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdu"
      f3d(n)%units = " "
      f3d(n)%data => hdu
      i_hdu = n ; n = n+1
!
      allocate(hdv(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('VN HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdv"
      f3d(n)%units = " "
      f3d(n)%data => hdv
      i_hdv = n ; n = n+1
!
      allocate(hdo2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O2 HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdo2"
      f3d(n)%units = " "
      f3d(n)%data => hdo2
      i_hdo2 = n ; n = n+1
!
      allocate(hdo1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('O1 HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdo1"
      f3d(n)%units = " "
      f3d(n)%data => hdo1
      i_hdo1 = n ; n = n+1
!
      allocate(hdhe(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('HE HORIZONTAL DIFFUSION')")
      f3d(n)%short_name = "hdhe"
      f3d(n)%units = " "
      f3d(n)%data => hdhe
      i_hdhe = n ; n = n+1
!
      allocate(ped(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('PEDERSEN CONDUCTIVITY')")
      f3d(n)%short_name = "ped"
      f3d(n)%units = " "
      f3d(n)%data => ped
      i_ped = n ; n = n+1
!
      allocate(hall(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('HALL CONDUCTIVITY')")
      f3d(n)%short_name = "hall"
      f3d(n)%units = " "
      f3d(n)%data => hall
      i_hall = n ; n = n+1
!
      allocate(lam1(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('LAMDA ION DRAG 1')")
      f3d(n)%short_name = "LAM1"
      f3d(n)%units = " "
      f3d(n)%data => lam1
      i_lam1 = n ; n = n+1
!
      allocate(zg(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('Geometric Height ZG')")
      f3d(n)%short_name = "ZG"
      f3d(n)%units = "cm"
      f3d(n)%data => zg
      i_zg = n ; n = n+1
!
      allocate(n2(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('Molecular Nitrogen N2')")
      f3d(n)%short_name = "N2"
      f3d(n)%units = "mmr"
      f3d(n)%data => n2
      i_n2 = n ; n = n+1
!
      allocate(wn(levd0:levd1,lond0:lond1,latd0:latd1),stat=istat)
      write(f3d(n)%long_name,"('Neutral Vertical Wind')")
      f3d(n)%short_name = "WN"
      f3d(n)%units = "cm"
      f3d(n)%data => wn
      i_wn = n ; n = n+1
!
! Check number of fields allocated:
      if (n-1 /= nf3d) then
        write(6,"(/,'>>> init_3d: wrong number of 3-d fields?',
     |    ' n-1=',i5,' nf3d=',i5)") n-1,nf3d
      endif
!
! f3d fields are diagnostics:
      f3d%prognostic = .false.
      f3d%magnetic = .false.
      f3d%task0_only = .false.
      f3d%vcoord = 'midpoints'
      do n=1,nf3d
        if (trim(f3d(n)%short_name)=='ZG') f3d(n)%vcoord='interfaces' 
      enddo
!
! Init:
      do n=1,nf3d
        f3d(n)%data = field_initval
      enddo
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 3-d fields:')") nf3d
        do n=1,nf3d
          call print_f3d(f3d(n),n)
        enddo
      endif
!
! Do other allocations (see allocdata.F):
      call allocdata
      end subroutine init_3d
!-----------------------------------------------------------------------
      subroutine init_fsech
!
! Initialize secondary history field structures. This does not include
! data, which is allocated and defined by sub addfld (addfld.F).
!
      use input_module,only: secflds
      use hist_module,only: nfsech
      implicit none
!
! Local:
      integer :: i,iprog,ier
!
! External:
      integer,external :: strloc
!
! Secondary history fields:
      fsechist(:)%long_name = ' '
      fsechist(:)%short_name= ' '
      fsechist(:)%units     = ' '
      fsechist(:)%prognostic  = .false.
      fsechist(:)%task0_only  = .false.
      fsechist(:)%mag = .false.
      fsechist(:)%geo = .false.

      do i=1,nfsech
        if (len_trim(secflds(i)) > 0) then
          iprog = strloc(f4d%short_name,nf4d,secflds(i))
!
! Is a prognostic: define fsech(i) from f4d(iprog):
!
          if (iprog > 0) then
            fsechist(i)%prognostic = .true.
            fsechist(i)%short_name = f4d(iprog)%short_name
            fsechist(i)%long_name  = f4d(iprog)%long_name
            fsechist(i)%units      = f4d(iprog)%units
!
! Is a diagnostic: define fsech(i)%name from input field name.
! Set units and long_name blank (will be optionally defined in
! user called sub addfld)
!
          else ! is diagnostic
            fsechist(i)%long_name = secflds(i)
            fsechist(i)%short_name = secflds(i)
            fsechist(i)%units = ' '
            fsechist(i)%long_name = ' '
          endif
        endif
      enddo ! i=1,nfsech
      end subroutine init_fsech
!-----------------------------------------------------------------------
      subroutine init_lbc(lon0,lon1,lat0,lat1)
      implicit none
!
! Arg:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Allocate t,u,v lbc (t,u,v at bottom interface level)
! These will be read from source history, then set in dt.F and duv.F.
!
! Subdomains:
      allocate(tlbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc')
      allocate(ulbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc')
      allocate(vlbc(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc')
      write(6,"('init_lbc: allocated subdomains tlbc, ulbc, vlbc')")
      
      allocate(tlbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_nm')
      allocate(ulbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_nm')
      allocate(vlbc_nm(lond0:lond1,latd0:latd1),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_nm')
      write(6,"('init_lbc_nm: allocated subdomains tlbc_nm, ulbc_nm,',
     |  ' vlbc_nm')")
!
! Global domain (subdomains will be gathered into these arrays
! for writing to the history, see sub mp_gather2root_lbc in mpi.F):
      allocate(tlbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_glb')
      allocate(ulbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_glb')
      allocate(vlbc_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_glb')
      write(6,"('init_lbc: allocated globals tlbc_glb, ulbc_glb,',
     |  ' vlbc_glb: nlonp4=',i3,' nlat=',i3)") nlonp4,nlat
     
      allocate(tlbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of tlbc_nm_glb')
      allocate(ulbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of ulbc_nm_glb')
      allocate(vlbc_nm_glb(nlonp4,nlat),stat=istat)
      if (istat /= 0) call shutdown('bad allocate of vlbc_nm_glb')
      write(6,"('init_lbc: allocated globals tlbc_nm_glb, ulbc_nm_glb,',
     |  ' vlbc_nm_glb: nlonp4=',i3,' nlat=',i3)") nlonp4,nlat
!
      end subroutine init_lbc
!-----------------------------------------------------------------------
      subroutine print_f4d(f,n)
      implicit none
      type(fields_4d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f4d
!-----------------------------------------------------------------------
      subroutine print_f3d(f,n)
      implicit none
      type(fields_3d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f3d
!-----------------------------------------------------------------------
      subroutine print_f2d(f,n)
      implicit none
      type(fields_2d),intent(in) :: f
      integer,intent(in) :: n
!
      write(6,"('Field ',i3,': ',a,' Short name: ',a,' Units: ',a)") 
     |  n,f%long_name(1:40),f%short_name(1:8),f%units(1:8)
      end subroutine print_f2d
!-----------------------------------------------------------------------
      end module fields_module
