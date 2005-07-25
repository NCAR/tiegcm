!
      module fields_module
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
     |  mxfsech,spval,nmlat,nmlonp1,nmagphrlon,nmagphrlat
      implicit none
      integer,parameter :: 
     |  nf4d          = 26, ! number of 4-d fields
     |  nf4d_hist     = 24, ! number of 4-d fields on primary histories
     |  nf3d          = 41  ! number of 3-d fields
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
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid
     |    magnetos,    ! true if field is on magnetospheric grid
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
     |  i_vc
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
     |  vc
!
! 3-d fields structure type:
! Data pointers will point to allocatable target arrays.
!
      type fields_3d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    magnetos,    ! true if field is on magnetospheric grid
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
     |  i_kldt   ,i_kldu    ,i_kldv    ,i_kldo2 ,i_kldo1 ,
     |  i_xnmbar ,i_xnmbari ,i_xnmbarm ,i_cp    ,i_kt    , 
     |  i_km     ,i_ui      ,i_vi      ,i_wi    ,i_vo2   , 
     |  i_vo1    ,i_vn2     ,i_sco2    ,i_sco1  ,i_scn2  ,
     |  i_xiop2p ,i_xiop2d  ,i_nplus   ,i_n2p   ,i_nop   ,
     |  i_lxx    ,i_lyy     ,i_lxy     ,i_lyx   ,i_qji_ti,
     |  i_qji_tn ,i_cool_implicit      ,i_cool_explicit  ,
     |  i_hdt    ,i_hdu     ,i_hdv     ,i_hdo2  ,i_hdo1  ,
     |  i_ped    ,i_hall    ,i_lam1
!
! Allocatable target arrays 3-d:
! Allocation will be:   tn(levd0:levd1,lond0 :lond1 ,latd0 :latd1 )
! which is the same as: tn(1    :nlev ,lon0-2:lon1+2,lat0-2:lat1+2)
!
      real,dimension(:,:,:),target,allocatable :: ! (k,i,j)
     |  kldt     ,kldu     ,kldv     ,kldo2    ,kldo1   ,
     |  xnmbar   ,xnmbari  ,xnmbarm  ,cp       ,kt      , 
     |  km       ,ui       ,vi       ,wi       ,vo2     , 
     |  vo1      ,vn2      ,sco2     ,sco1     ,scn2    ,
     |  xiop2p   ,xiop2d   ,nplus    ,n2p      ,nop     ,
     |  lxx      ,lyy      ,lxy      ,lyx      ,qji_ti  ,
     |  qji_tn   ,cool_implicit      ,cool_explicit     ,
     |  hdt      ,hdu      ,hdv      ,hdo2     ,hdo1    ,
     |  ped      ,hall     ,lam1
!
! 2-d field type (used for 2d secondary history fields):
      type fields_2d
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        logical                      :: 
     |    magnetic,    ! true if field is on magnetic 3d grid 
     |    magnetos,    ! true if field is on magnetospheric grid
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
     |  dynpot(nlonp1,0:nlatp1,nlevp1),  ! 3d electric potential geographic
     |  phim3d(nmlonp1,nmlat,-2:nlevp1)  ! 3d electric potential magnetic
!
! Secondary history fields (runtime user defined 3-d fields):
!
      type(fields_3d) :: fsech(mxfsech)         ! geographic
      type(fields_3d) :: fsechmag(mxfsech)      ! geomagnetic
      type(fields_2d) :: fsech2d(mxfsech)       ! 2-d geographic 
      type(fields_2d) :: fsechmag2d(mxfsech)    ! 2-d magnetic
      type(fields_2d) :: fsechmagphr2d(mxfsech) ! 2-d magnetosphere
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
! Full 3d grid with all primary history fields for writing to netcdf
! history files. This will be allocated only on the root task 
! (see allocdata):
      real,allocatable :: foutput(:,:,:,:) ! (nlevp1,nlonp4,nlat,nf4d_hist)
!
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
      integer :: n,istat
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
      f4d(n)%data => tn
      i_tn = n ; n = n+1
!
      allocate(un(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND')")
      f4d(n)%short_name = "UN"
      f4d(n)%units = "CM/S"
      f4d(n)%data => un
      i_un = n ; n = n+1
!
      allocate(vn(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND')")
      f4d(n)%short_name = "VN"
      f4d(n)%units = "CM/S"
      f4d(n)%data => vn
      i_vn = n ; n = n+1
!
      allocate(o2(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN')")
      f4d(n)%short_name = ("O2")
      f4d(n)%units = "MMR"
      f4d(n)%data => o2
      i_o2 = n ; n = n+1
!
      allocate(o1(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN')")
      f4d(n)%short_name = ("O1")
      f4d(n)%units = "MMR"
      f4d(n)%data => o1
      i_o1 = n ; n = n+1
!
      allocate(n4s(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S')")
      f4d(n)%short_name = ("N4S")
      f4d(n)%units = "MMR"
      f4d(n)%data => n4s
      i_n4s = n ; n = n+1
!
      allocate(no(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NITRIC OXIDE')")
      f4d(n)%short_name = ("NO")
      f4d(n)%units = "MMR"
      f4d(n)%data => no
      i_no = n ; n = n+1
!
      allocate(op(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O+ ION')")
      f4d(n)%short_name = ("OP")
      f4d(n)%units = "CM^3"
      f4d(n)%data => op
      i_op = n ; n = n+1
!
      allocate(n2d(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N2D')")
      f4d(n)%short_name = ("N2D")
      f4d(n)%units = "MMR"
      f4d(n)%data => n2d
      i_n2d = n ; n = n+1
!
      allocate(ti(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ION TEMPERATURE')")
      f4d(n)%short_name = ("TI")
      f4d(n)%units = "DEG K"
      f4d(n)%data => ti
      i_ti = n ; n = n+1
!
      allocate(te(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRON TEMPERATURE')")
      f4d(n)%short_name = ("TE")
      f4d(n)%units = "DEG K"
      f4d(n)%data => te
      i_te = n ; n = n+1
!
      allocate(ne(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRON DENSITY')")
      f4d(n)%short_name = ("NE")
      f4d(n)%units = "CM^3"
      f4d(n)%data => ne
      i_ne = n ; n = n+1
!
      allocate(o2p(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('O2+ ION')")
      f4d(n)%short_name = ("O2P")
      f4d(n)%units = "CM^3"
      f4d(n)%data => o2p
      i_o2p = n ; n = n+1
!
      allocate(w(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('VERTICAL VELOCITY (PLUS UP)')")
      f4d(n)%short_name = ("W")
      f4d(n)%units = "CM/S"
      f4d(n)%data => w
      i_w = n ; n = n+1
!
      allocate(z(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('GEOPOTENTIAL HEIGHT')")
      f4d(n)%short_name = ("Z")
      f4d(n)%units = "CM"
      f4d(n)%data => z
      i_z = n ; n = n+1
!
      allocate(poten(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ELECTRIC POTENTIAL')")
      f4d(n)%short_name = ("POTEN")
      f4d(n)%units = "VOLTS"
      f4d(n)%data => poten
      i_poten = n ; n = n+1
!
      allocate(tn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL TEMPERATURE (TIME N-1)')")
      f4d(n)%short_name = ("TN_NM")
      f4d(n)%units = "DEG K"
      f4d(n)%data => tn_nm
      i_tn_nm = n ; n = n+1
!
      allocate(un_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL ZONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("UN_NM")
      f4d(n)%units = "CM/S"
      f4d(n)%data => un_nm
      i_un_nm = n ; n = n+1
!
      allocate(vn_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NEUTRAL MERIDIONAL WIND (TIME N-1)')")
      f4d(n)%short_name = ("VN_NM")
      f4d(n)%units = "CM/S"
      f4d(n)%data => vn_nm
      i_vn_nm = n ; n = n+1
!
      allocate(o2_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MOLECULAR OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O2_NM")
      f4d(n)%units = "MMR"
      f4d(n)%data => o2_nm
      i_o2_nm = n ; n = n+1
!
      allocate(o1_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('ATOMIC OXYGEN (TIME N-1)')")
      f4d(n)%short_name = ("O1_NM")
      f4d(n)%units = "MMR"
      f4d(n)%data => o1_nm
      i_o1_nm = n ; n = n+1
!
      allocate(n4s_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('N4S (TIME N-1)')")
      f4d(n)%short_name = ("N4S_NM")
      f4d(n)%units = "MMR"
      f4d(n)%data => n4s_nm
      i_n4s_nm = n ; n = n+1
!
      allocate(no_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('NO (TIME N-1)')")
      f4d(n)%short_name = ("NO_NM")
      f4d(n)%units = "MMR"
      f4d(n)%data => no_nm
      i_no_nm = n ; n = n+1
!
      allocate(op_nm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('OP (TIME N-1)')")
      f4d(n)%short_name = ("OP_NM")
      f4d(n)%units = "CM^3"
      f4d(n)%data => op_nm
      i_op_nm = n ; n = n+1
!
      allocate(barm(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('MEAN MOLECULAR WEIGHT')")
      f4d(n)%short_name = ("BARM")
      f4d(n)%units = ' '
      f4d(n)%data => barm
      i_barm = n ; n = n+1
!
      allocate(vc(levd0:levd1,lond0:lond1,latd0:latd1,2),stat=istat)
      write(f4d(n)%long_name,"('COS(PHI)*VN')")
      f4d(n)%short_name = ("VC")
      f4d(n)%units = ' '
      f4d(n)%data => vc
      i_vc = n ; n = n+1
!
! itp,itc are indices to rightmost dimension of field data, for
! previous and current time steps.
      itp = 1
      itc = 2
!
! Set polesign for crossing the poles (see mk_polelats and mp_bndlats)
      f4d%polesign = 1.
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
      f4d%magnetos = .false.
!
! f4d fields are "prognostic":
      f4d%prognostic = .true.
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 4-d fields:')") nf4d
        do n=1,nf4d
          call print_f4d(f4d(n),n)
        enddo
      endif
!
! Report to stdout:
      if (iprint > 0) then
        write(6,"(/,'There are ',i3,' 4-d fields:')") nf4d
        do n=1,nf4d
          call print_f4d(f4d(n),n)
        enddo
      endif
!
! Do other allocations:
!     call allocdata
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
! f3d fields are diagnostics:
      f3d%prognostic = .false.
      f3d%magnetic = .false.
      f3d%magnetos = .false.
      f3d%task0_only = .false.
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
      subroutine set_fsech
      use input_module,only: secflds,secfmag,secfgeo2d,secfmag2d,
     |   secfmagphr
      use hist_module,only: nfsech_geo,nfsech_mag,nfsech_geo2d,
     |   nfsech_mag2d,nfsech_magphr
!
! Initialize secondary history fields (to be defined at runtime by 
! user-callable subroutine addfsech):
!
! Local:
      integer :: i,iprog,idiag,ier
!
! External:
      integer,external :: strloc
!
! Secondary history fields on geographic grid:
      fsech(:)%long_name = ' '
      fsech(:)%short_name= ' '
      fsech(:)%units     = ' '
      fsech(:)%magnetic  = .false.
      fsech(:)%magnetos  = .false.
!
      do i=1,nfsech_geo
        if (len_trim(secflds(i)) > 0) then
          iprog = strloc(f4d%short_name,nf4d,secflds(i))
!
! Is a prognostic: define fsech(i) from f4d(iprog):
!
          if (iprog > 0) then
            fsech(i)%prognostic = .true.
            fsech(i)%short_name = f4d(iprog)%short_name
            fsech(i)%long_name  = f4d(iprog)%long_name
            fsech(i)%units      = f4d(iprog)%units
!
! Is a diagnostic: define fsech(i)%name from input field name.
! Set units and long_name blank (will be optionally defined in
! user called sub addfsech)
!
          else ! is diagnostic
            fsech(i)%long_name = secflds(i)
            fsech(i)%short_name = secflds(i)
            fsech(i)%units = ' '
            fsech(i)%long_name = ' '
          endif
!
! Allocate pointer to data for 3d diagnostic field and initialize
! the data to spval. The field should be defined later by user-called
! sub addfsech. 
!
          allocate(fsech(i)%data(nlevp1,nlonp4,nlat),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> WARNING set_fsech:',
     |        ' error allocating fsech(i)%data: i=',i3)") i
          else
            write(6,"('Allocated (nlevp1=',i3,
     |        ',nlonp4=',i3,',nlat=',i3,') for geographic field ',a)") 
     |        nlevp1,nlonp4,nlat,trim(fsech(i)%short_name)
          endif
          fsech(i)%data = spval ! array op
        endif
      enddo 
!
! 3d secondary history fields on magnetic grid:
      fsechmag(:)%long_name = ' '
      fsechmag(:)%short_name= ' '
      fsechmag(:)%units     = ' '
      fsechmag(:)%magnetic  = .true.
      fsechmag(:)%magnetos  = .false.
      fsechmag(:)%prognostic = .true.
      fsechmag(:)%task0_only = .false.
      do i=1,nfsech_mag
        if (len_trim(secfmag(i)) > 0) then
          fsechmag(i)%long_name = secfmag(i)
          fsechmag(i)%short_name = secfmag(i)
          allocate(fsechmag(i)%data(nlevp1+3,nmlonp1,nmlat),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> WARNING set_fsech:',
     |        ' error allocating fsechmag(i)%data: i=',i3)") i
          else
            write(6,"('Allocated sechist data(nlevp1+3=',i3,
     |        ' nmlonp1=',i3,',nmlat=',i3,') for mag field ',a)") 
     |        nlevp1+3,nmlonp1,nmlat,trim(fsechmag(i)%short_name)
          endif
          fsechmag(i)%data = spval ! array op
        endif
      enddo
!
! 2d secondary history fields geographic fields:
      fsech2d(:)%long_name = ' '
      fsech2d(:)%short_name= ' '
      fsech2d(:)%units     = ' '
      fsech2d(:)%magnetic  = .false.  ! will be determined from addfsech
      fsech2d(:)%magnetos  = .false.
      fsech2d(:)%prognostic = .false.
      fsech2d(:)%task0_only = .false.
      do i=1,nfsech_geo2d
        if (len_trim(secfgeo2d(i)) > 0) then
          fsech2d(i)%long_name = secfgeo2d(i)
          fsech2d(i)%short_name = secfgeo2d(i)
          allocate(fsech2d(i)%data(nlonp4,nlat),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> WARNING set_fsech:',
     |        ' error allocating fsech2d(i)%data: i=',i3)") i
          else
            write(6,"('Allocated sechist data(nlonp4=',i3,
     |        ' nlat=',i3,') for geo2d field ',a)") 
     |        nlonp4,nlat,trim(fsech2d(i)%short_name)
          endif
          fsech2d(i)%data = spval ! array op
        endif
      enddo
!
! 2d secondary history magnetic fields:
      fsechmag2d(:)%long_name = ' '
      fsechmag2d(:)%short_name= ' '
      fsechmag2d(:)%units     = ' '
      fsechmag2d(:)%magnetic  = .true.
      fsechmag2d(:)%magnetos  = .false.
      fsechmag2d(:)%prognostic = .false.
      fsechmag2d(:)%task0_only = .false.
      do i=1,nfsech_mag2d
        if (len_trim(secfmag2d(i)) > 0) then
          fsechmag2d(i)%long_name  = secfmag2d(i)
          fsechmag2d(i)%short_name = secfmag2d(i)
          allocate(fsechmag2d(i)%data(nmlonp1,nmlat),stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> WARNING set_fsech:',
     |        ' error allocating fsechmag2d(i)%data: i=',i3)") i
          else
            write(6,"('Allocated sechist data(nmlonp1=',i3,
     |        ' nmlat=',i3,') for mag2d field ',a)") 
     |        nmlonp1,nmlat,trim(fsechmag2d(i)%short_name)
          endif
          fsechmag2d(i)%data = spval ! array op
        endif
      enddo
!
! 2d secondary history magnetospheric fields:
      fsechmagphr2d(:)%long_name = ' '
      fsechmagphr2d(:)%short_name= ' '
      fsechmagphr2d(:)%units     = ' '
      fsechmagphr2d(:)%magnetic  = .false.
      fsechmagphr2d(:)%magnetos  = .true.
      fsechmagphr2d(:)%prognostic= .false.
      fsechmagphr2d(:)%task0_only= .false.
      do i=1,nfsech_magphr
        if (len_trim(secfmagphr(i)) > 0) then
          fsechmagphr2d(i)%long_name  = secfmagphr(i)
          fsechmagphr2d(i)%short_name = secfmagphr(i)
          allocate(fsechmagphr2d(i)%data(nmagphrlon,nmagphrlat),
     |       stat=ier)
          if (ier /= 0) then
            write(6,"(/,'>>> WARNING set_fsech:',
     |        ' error allocating fsechmagphr2d(i)%data: i=',i3)") i
          else
            write(6,"('Allocated sechist data(nmagphrlon=',i3,
     |        ' nmagphrlat=',i3,') for magphr field ',a)") 
     |        nmagphrlon,nmagphrlat,trim(fsechmagphr2d(i)%short_name)
          endif
          fsechmagphr2d(i)%data = spval ! array op
        endif
      enddo
      end subroutine set_fsech
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
