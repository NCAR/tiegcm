      module diags_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Make "sanctioned" diagnostics available on secondary history files.
! To request field(s), the user adds Shortname(s) to the namelist SECFLDS.
!
! A table of available diagnostic fields, including names, levels, units, 
! etc., is printed to stdout by sub init_diags.
!
      use mpi_module,only: mytid
      use addfld_module,only: addfld
      use params_module,only: dlev
      implicit none
      integer,parameter :: 
     |  longname_len  = 80, ! length of field long name 
     |  shortname_len = 16, ! length of field short name 
     |  units_len     = 16  ! length of field units attribute
!
! ndiag: Number of available diagnostic fields. This must be incremented when 
! adding new fields, and type diags(ndiag) must be defined by sub init_diags 
! below.
!
      integer,parameter :: ndiag = 43 ! increment this to add new diags
      type diag_type
        character(len=longname_len)  :: long_name
        character(len=shortname_len) :: short_name
        character(len=units_len)     :: units
        character(len=shortname_len) :: levels   ! "lev", "ilev", or "none"
        character(len=shortname_len) :: caller   ! File or sub that calls mkdiag_xxx
      end type diag_type
      type(diag_type) :: diags(ndiag)

      contains
!-----------------------------------------------------------------------
      subroutine init_diags(iprint)
      implicit none
!
! Args:
      integer,intent(in) :: iprint ! if > 0, print diags list to stdout
!
! Init to empty strings:
!
      integer :: n
      do n=1,ndiag
        diags(n)%long_name  = ' '
        diags(n)%short_name = ' '
        diags(n)%units      = ' '
        diags(n)%levels     = ' '
      enddo
!
! Fields order is approximately (not necesarilly) alphabetical by short_name:
! Table of these values are written to stdout, see below.
!
      n = 1
      diags(n)%short_name = 'CO2_COOL'
      diags(n)%long_name  = 'CO2 Cooling'
      diags(n)%units      = 'erg/g/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'newton.F'
!
      n = n+1
      diags(n)%short_name = 'NO_COOL'
      diags(n)%long_name  = 'NO Cooling'
      diags(n)%units      = 'erg/g/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'newton.F'
!
      n = n+1
      diags(n)%short_name = 'DEN'
      diags(n)%long_name  = 'Total Density'
      diags(n)%units      = 'g/cm3'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'dt.F'
!
      n = n+1
      diags(n)%short_name = 'HEATING'
      diags(n)%long_name  = 'Total Heating'
      diags(n)%units      = 'erg/g/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'dt.F'
!
      n = n+1
      diags(n)%short_name = 'HMF2'
      diags(n)%long_name  = 'HMF2: Height of the F2 Layer' 
      diags(n)%units      = 'km'
      diags(n)%levels     = 'none' ! hmf2 is 2d lon x lat
      diags(n)%caller     = 'elden.F'
!
      n = n+1
      diags(n)%short_name = 'NMF2'
      diags(n)%long_name  = 'NMF2: Peak Density of the F2 Layer' 
      diags(n)%units      = '1/cm3'
      diags(n)%levels     = 'none' ! nmf2 is 2d lon x lat
      diags(n)%caller     = 'elden.F'
!
      n = n+1
      diags(n)%short_name = 'FOF2'
      diags(n)%long_name  = 'FOF2: Critical Frequency of the F2 Layer' 
      diags(n)%units      = 'MHz'
      diags(n)%levels     = 'none' ! fof2 is 2d lon x lat
      diags(n)%caller     = 'elden.F'
!
      n = n+1
      diags(n)%short_name = 'JE13D'
      diags(n)%long_name  = 'JE13D: Eastward current density (3d)'
      diags(n)%units      = 'A/m2'
      diags(n)%levels     = 'mlev' 
      diags(n)%caller     = 'current.F90'
!
      n = n+1
      diags(n)%short_name = 'JE23D'
      diags(n)%long_name  = 'JE23D: Downward current density (3d)'
      diags(n)%units      = 'A/m2'
      diags(n)%levels     = 'mlev' 
      diags(n)%caller     = 'current.F90'
!     
      n = n+1
      diags(n)%short_name = 'JQR'
      diags(n)%long_name  = 'JQR: Upward current density (2d)'
      diags(n)%units      = 'A/m2'
      diags(n)%levels     = 'none' 
      diags(n)%caller     = 'current.F90'
!     
      n = n+1
      diags(n)%short_name = 'KQLAM'
      diags(n)%long_name  = 
     |  'KQLAM: Height-integrated current density (+north)'
      diags(n)%units      = 'A/m'
      diags(n)%levels     = 'none' 
      diags(n)%caller     = 'current.F90'
!     
      n = n+1
      diags(n)%short_name = 'KQPHI'
      diags(n)%long_name  = 
     |  'KQPHI: Height-integrated current density (+east)'
      diags(n)%units      = 'A/m'
      diags(n)%levels     = 'none' 
      diags(n)%caller     = 'current.F90'
!
      n = n+1
      diags(n)%short_name = 'LAMDA_HAL'
      diags(n)%long_name  = 'LAMDA_HAL: Hall Ion Drag Coefficient'
      diags(n)%units      = '1/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'lamdas.F'
!
      n = n+1
      diags(n)%short_name = 'LAMDA_PED'
      diags(n)%long_name  = 'LAMDA_PED: Pedersen Ion Drag Coefficient'
      diags(n)%units      = '1/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'lamdas.F'
!
      n = n+1
      diags(n)%short_name = 'MU_M'
      diags(n)%long_name  = 'MU_M: Molecular Viscosity Coefficient'
      diags(n)%units      = 'g/cm/s'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'cpktkm.F'
!
      n = n+1
      diags(n)%short_name = 'QJOULE'
      diags(n)%long_name  = 'QJOULE: Joule Heating'
      diags(n)%units      = 'erg/g/s'  ! model description has erg/K/s
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'qjoule.F'
!
      n = n+1
      diags(n)%short_name = 'SCHT'
      diags(n)%long_name  = 'SCHT: Pressure Scale Height'
      diags(n)%units      = 'km'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'addiag.F'
!
      n = n+1
      diags(n)%short_name = 'SIGMA_HAL'
      diags(n)%long_name  = 'SIGMA_HAL: Hall Conductivity'
      diags(n)%units      = 'S/m'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'lamdas.F'
!
      n = n+1
      diags(n)%short_name = 'SIGMA_PED'
      diags(n)%long_name  = 'SIGMA_PED: Pedersen Conductivity'
      diags(n)%units      = 'S/m'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'lamdas.F'
!
      n = n+1
      diags(n)%short_name = 'TEC'
      diags(n)%long_name  = 'TEC: Total Electron Content'
      diags(n)%units      = '1/cm2'
      diags(n)%levels     = 'none' ! 2d lon x lat
      diags(n)%caller     = 'elden.F'
!
      n = n+1
      diags(n)%short_name = 'UI_ExB'
      diags(n)%long_name  = 'UI: Zonal ExB Velocity'
      diags(n)%units      = 'cm/s'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'ionvel.F'
!
      n = n+1
      diags(n)%short_name = 'VI_ExB'
      diags(n)%long_name  = 'VI: Meridional ExB Velocity'
      diags(n)%units      = 'cm/s'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'ionvel.F'
!
      n = n+1
      diags(n)%short_name = 'WI_ExB'
      diags(n)%long_name  = 'WI: Vertical ExB Velocity'
      diags(n)%units      = 'cm/s'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'ionvel.F'
!
      n = n+1
      diags(n)%short_name = 'WN'
      diags(n)%long_name  = 'WN: Neutral Vertical Wind (plus up)' 
      diags(n)%units      = 'cm/s'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'swdot.F'
!
      n = n+1
      diags(n)%short_name = 'O_N2'
      diags(n)%long_name  = 'O/N2 RATIO'
      diags(n)%units      = 'none'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'comp.F'
!
      n = n+1
      diags(n)%short_name = 'QJOULE_INTEG'
      diags(n)%long_name  = 'Height-integrated Joule Heating'
      diags(n)%units      = 'erg/cm2/s'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'qjoule.F'
!
      n = n+1
      diags(n)%short_name = 'BX'
      diags(n)%long_name  = 
     |  'BX/BMAG: Normalized eastward component of magnetic field'
      diags(n)%units      = ' '
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'oplus.F'
!
      n = n+1
      diags(n)%short_name = 'BY'
      diags(n)%long_name  = 
     |  'BY/BMAG: Normalized northward component of magnetic field'
      diags(n)%units      = ' '
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'oplus.F'
!
      n = n+1
      diags(n)%short_name = 'BZ'
      diags(n)%long_name  = 
     |  'BZ/BMAG: Normalized upward component of magnetic field'
      diags(n)%units      = ' '
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'oplus.F'
!
      n = n+1
      diags(n)%short_name = 'BMAG'
      diags(n)%long_name  = 'BMAG: Magnetic field magnitude'
      diags(n)%units      = 'Gauss'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'oplus.F'
!
      n = n+1
      diags(n)%short_name = 'EX'
      diags(n)%long_name  = 
     |  'EX: Zonal component of electric field'
      diags(n)%units      = 'V/m'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'EY'
      diags(n)%long_name  = 
     |  'EY: Meridional component of electric field'
      diags(n)%units      = 'V/m'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'EZ'
      diags(n)%long_name  = 
     |  'EZ: Vertical component of electric field'
      diags(n)%units      = 'V/m'
      diags(n)%levels     = 'ilev'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'ED1'
      diags(n)%long_name  = 
     |  'ED1: Magnetic eastward component of electric field'
      diags(n)%units      = 'V/m'
      diags(n)%levels     = 'imlev'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'ED2'
      diags(n)%long_name  = 
     |'ED2: Magnetic downward (equatorward) component of electric field'
      diags(n)%units      = 'V/m'
      diags(n)%levels     = 'imlev'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'PHIM2D'
      diags(n)%long_name  = 
     |  'PHIM2D: 2d Electric Potential on magnetic grid'
      diags(n)%units      = 'V'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'pdynamo.F'
!
      n = n+1
      diags(n)%short_name = 'N2'
      diags(n)%long_name  = 'N2: Molecular Nitrogen'
      diags(n)%units      = 'mmr'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'comp.F'
!
      n = n+1
      diags(n)%short_name = 'ZGMID'
      diags(n)%long_name  = 'ZGMID: Geometric Height at midpoints'
      diags(n)%units      = 'cm'
      diags(n)%levels     = 'lev'
      diags(n)%caller     = 'addiag.F'
!
      n = n+1
      diags(n)%short_name = 'CUSP'
      diags(n)%long_name  = 'CUSP Low Energy Electron Flux'
      diags(n)%units      = 'erg/cm2/s'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'dynamics.F'
!
      n = n+1
      diags(n)%short_name = 'DRIZZLE'
      diags(n)%long_name  = 'DRIZZLE Low Energy Electron Flux'
      diags(n)%units      = 'erg/cm2/s'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'dynamics.F'
!
      n = n+1
      diags(n)%short_name = 'ALFA'
      diags(n)%long_name  = 'Aurora Characteristic Energy'
      diags(n)%units      = 'keV'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'dynamics.F'
!
      n = n+1
      diags(n)%short_name = 'NFLUX'
      diags(n)%long_name  = 'Aurora Number Flux'
      diags(n)%units      = '#/cm2/s'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'dynamics.F'
!
      n = n+1
      diags(n)%short_name = 'EFLUX'
      diags(n)%long_name  = 'Aurora Energy Flux'
      diags(n)%units      = 'erg/cm2/s'
      diags(n)%levels     = 'none'
      diags(n)%caller     = 'dynamics.F'

      if (n /= ndiag) then
        write(6,"('>>> init_diags: n=',i4,' ndiag=',i4)")
     |    n,ndiag
        call shutdown('init_diags ndiag')
      endif
!
! Print table to stdout:
      if (iprint > 0) then
        write(6,"(/,96('-'))")
        write(6,"('Table of Available Diagnostic Fields:')")
        write(6,"('Shortnames may be added to namelist SECFLDS',/)")
        write(6,"('Field   Shortname       Units           Levels  ',
     |    'Caller          Longname')")
        write(6,"(48('- '))")
        do n=1,ndiag
          write(6,"(i4,4x,a,a,a,a,a)") n,diags(n)%short_name,
     |      diags(n)%units,diags(n)%levels(1:8),diags(n)%caller,
     |      trim(diags(n)%long_name)
        enddo
        write(6,"(96('-'),/)")
      endif
      end subroutine init_diags
!-----------------------------------------------------------------------
      integer function checkf(name)
!
! Return index to name in diags(:).
!
      use input_module,only: secflds,mxfsech
!
! Args:
      character(len=*),intent(in) :: name
!
! External:
      integer,external :: strloc ! util.F
!
! Check for field in namelist SECFLDS -- if not found, return silently
      checkf = strloc(secflds,mxfsech,name)
      if (checkf==0) return
!
! Get index to diags -- if not found, return zero with error msg:
      checkf = strloc(diags(:)%short_name,ndiag,name)
      if (checkf==0) then
        write(6,"('>>> Cannot find field ',a,' in list of ',
     |    'available diagnostics -- returning')") name
      endif
      return
      end function checkf
!-----------------------------------------------------------------------
      subroutine mkdiag_DEN(name,rho,lev0,lev1,lon0,lon1,lat)
!
! Save diagnostic DEN (Total Density). This routine is called from 
! dt.F, and simply passes rho to addfld at the current latitude.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: rho(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save to secondary history. Note rho is calculated on interface
! levels, hence 'ilev' in the addfld call:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,rho,diags(ix)%levels,lev0,lev1,
     |  'lon',lon0,lon1,lat)
      end subroutine mkdiag_DEN
!-----------------------------------------------------------------------
      subroutine mkdiag_QJOULE(name,qji_tn,lev0,lev1,lon0,lon1,lat)
!
! From qjoule.F, levels at midpoints.
!
! Args
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: qji_tn(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary history:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,qji_tn(lev0:lev1-1,lon0:lon1),diags(ix)%levels,
     |  lev0,lev1-1,'lon',lon0,lon1,lat)

      end subroutine mkdiag_QJOULE
!-----------------------------------------------------------------------
      subroutine mkdiag_TEC(name,tec,z,electrons,lev0,lev1,lon0,lon1,
     |  lat)
!
! Calculate 2d global height-integrated diagnostic TEC (Total Electron Content)
! This is called from sub elden, elden.F at latitude lat.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: z,electrons
      real,intent(out) :: tec(lon0:lon1)
!
! Local:
      integer :: i,k
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Integrate electron content in height at current latitude:
      tec(:) = 0.
      do i=lon0,lon1
        do k=lev0,lev1-1
          tec(i) = tec(i)+(z(k+1,i)-z(k,i))*electrons(k,i)
        enddo
      enddo
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,tec(:),'lon',lon0,lon1,'lat',lat,lat,0)
      end subroutine mkdiag_TEC
!-----------------------------------------------------------------------
      subroutine mkdiag_UI(name,ui,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic UI_ExB. This is called from ionvel.F. 
! This routine simply passes ui from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: ui(lev0:lev1,lon0:lon1,lat0:lat1)
!
! Local:
      integer :: lat
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories (Convert back to cm/s):
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, ui(:,:,lat),diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_UI
!-----------------------------------------------------------------------
      subroutine mkdiag_VI(name,vi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic VI_ExB. This is called from ionvel.F. 
! This routine simply passes vi from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: vi(lev0:lev1,lon0:lon1,lat0:lat1)
!
! Local:
      integer :: lat
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, vi(:,:,lat),diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_VI
!-----------------------------------------------------------------------
      subroutine mkdiag_WI(name,wi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save 3d diagnostic WI_ExB. This is called from ionvel.F. 
! This routine simply passes wi from ionvel to addfld.
! Levels at interfaces. Incoming units are cm/s, but is put
! on the history in m/s.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: wi(lev0:lev1,lon0:lon1,lat0:lat1)
!
! Local:
      integer :: lat
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      do lat=lat0,lat1
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units, wi(:,:,lat),diags(ix)%levels,
     |    lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      end subroutine mkdiag_WI
!-----------------------------------------------------------------------
      subroutine mkdiag_HEAT(name,total_heat,lev0,lev1,lon0,lon1,lat)
!
! Total heating from dt.F. Levels on midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: total_heat(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,total_heat(lev0:lev1-1,:),diags(ix)%levels,
     |  lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_HEAT
!-----------------------------------------------------------------------
      subroutine mkdiag_CO2COOL(name,co2_cool,lev0,lev1,lon0,lon1,lat)
!
! CO2 cooling from newton.F, levels="lev" (midpoints)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: co2_cool(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,co2_cool,diags(ix)%levels,lev0,lev1,
     |  'lon',lon0,lon1,lat)
      end subroutine mkdiag_CO2COOL
!-----------------------------------------------------------------------
      subroutine mkdiag_NOCOOL(name,no_cool,lev0,lev1,lon0,lon1,lat)
!
! NO cooling from newton.F, levels="lev" (midpoints)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: no_cool(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,no_cool,diags(ix)%levels,lev0,lev1,
     |  'lon',lon0,lon1,lat)
      end subroutine mkdiag_NOCOOL
!-----------------------------------------------------------------------
      subroutine mkdiag_SIGMAPED(name,sigmaped,lev0,lev1,lon0,lon1,lat)
!
! From lamdas.F. Levels at midpoints (altho lamdas in lamdas.F are
! at interfaces)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: sigmaped(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,sigmaped(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_SIGMAPED
!-----------------------------------------------------------------------
      subroutine mkdiag_SIGMAHAL(name,sigmahal,lev0,lev1,lon0,lon1,lat)
!
! From lamdas.F. Levels at midpoints (altho lamdas themselves are at
! interfaces)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: sigmahal(lev0:lev1,lon0-2:lon1+2)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,sigmahal(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_SIGMAHAL
!-----------------------------------------------------------------------
      subroutine mkdiag_LAMDAPED(name,lamdaped,lev0,lev1,lon0,lon1,lat)
!
! Save Pedersen ion drag coefficient (on interfaces, called from lamdas.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: lamdaped(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,lamdaped(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_LAMDAPED
!-----------------------------------------------------------------------
      subroutine mkdiag_LAMDAHAL(name,lamdahal,lev0,lev1,lon0,lon1,lat)
!
! Save Hall ion drag coefficient (on interfaces, called from lamdas.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in) :: lamdahal(lev0:lev1,lon0:lon1)
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Save on secondary histories:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,lamdahal(lev0:lev1-1,lon0:lon1),
     |  diags(ix)%levels,lev0,lev1-1,'lon',lon0,lon1,lat)
      end subroutine mkdiag_LAMDAHAL
!-----------------------------------------------------------------------
      subroutine mkdiag_EXYZ(name,exyz,lev0,lev1,lon0,lon1,lat)
!
! 3d electric field on geographic grid (called by efield.F)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in) :: exyz
!
! Local:
      integer,save :: iex=-1,iey=-1,iez=-1

      if (trim(name)=='EX'.and.iex==-1) iex = checkf(name)
      if (trim(name)=='EY'.and.iey==-1) iey = checkf(name)
      if (trim(name)=='EZ'.and.iez==-1) iez = checkf(name)
      if (iex<=0.and.iey<=0.and.iez<=0) return
!
! Save 2d field (3d on history) at current lat on secondary histories:
! Pass exyz*100. to convert from V/cm to V/m
!
      if (trim(name)=='EX'.and.iex > 0) then
        call addfld(diags(iex)%short_name,diags(iex)%long_name,
     |    diags(iex)%units,exyz*100.,diags(iex)%levels,lev0,lev1,'lon',
     |    lon0,lon1,lat)
      elseif (trim(name)=='EY'.and.iey > 0) then
        call addfld(diags(iey)%short_name,diags(iey)%long_name,
     |    diags(iey)%units,exyz*100.,diags(iey)%levels,lev0,lev1,'lon',
     |    lon0,lon1,lat)
      elseif (trim(name)=='EZ'.and.iez > 0) then
        call addfld(diags(iez)%short_name,diags(iez)%long_name,
     |    diags(iez)%units,exyz*100.,diags(iez)%levels,lev0,lev1,'lon',
     |    lon0,lon1,lat)
      endif
      end subroutine mkdiag_EXYZ
!-----------------------------------------------------------------------
      subroutine mkdiag_BXYZ(name,bxyz,lon0,lon1,lat)
!
! Called by oplus.F
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat
      real,dimension(lon0:lon1),intent(in) :: bxyz
!
! Local:
      integer,save :: ibx=-1,iby=-1,ibz=-1

      if (trim(name)=='BX'.and.ibx==-1) ibx = checkf(name)
      if (trim(name)=='BY'.and.iby==-1) iby = checkf(name)
      if (trim(name)=='BZ'.and.ibz==-1) ibz = checkf(name)
      if (ibx<=0.and.iby<=0.and.ibz<=0) return
!
! Save 2d field on secondary histories:
      if (trim(name)=='BX'.and.ibx > 0) then
        call addfld(diags(ibx)%short_name,diags(ibx)%long_name,
     |    diags(ibx)%units,bxyz,'lon',lon0,lon1,'lat',lat,lat,0)
      elseif (trim(name)=='BY'.and.iby > 0) then
        call addfld(diags(iby)%short_name,diags(iby)%long_name,
     |    diags(iby)%units,bxyz,'lon',lon0,lon1,'lat',lat,lat,0)
      elseif (trim(name)=='BZ'.and.ibz > 0) then
        call addfld(diags(ibz)%short_name,diags(ibz)%long_name,
     |    diags(ibz)%units,bxyz,'lon',lon0,lon1,'lat',lat,lat,0)
      endif
      end subroutine mkdiag_BXYZ
!-----------------------------------------------------------------------
      subroutine mkdiag_BMAG(name,bmag,lon0,lon1,lat)
!
! Called by oplus.F
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat
      real,dimension(lon0:lon1),intent(in) :: bmag
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,bmag,'lon',lon0,lon1,'lat',lat,lat,0)

      end subroutine mkdiag_BMAG
!-----------------------------------------------------------------------
      subroutine mkdiag_ED12(name,ed12,mlon0,mlon1,mlev0,mlev1,lat)
!
! Two 3d components of electric field on the geomagnetic grid: ED1,ED2 
! (called by dynamo.F with global magnetic arrays, mlev0:mlev1
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: mlon0,mlon1,mlev0,mlev1,lat
      real,dimension(mlon0:mlon1,mlev0:mlev1),intent(in) :: ed12
!
! Local:
      integer,save :: ied1=-1,ied2=-1

      if (trim(name)=='ED1'.and.ied1==-1) ied1 = checkf(name)
      if (trim(name)=='ED2'.and.ied2==-1) ied2 = checkf(name)
      if (ied1<=0.and.ied2<=0) return
!
! Save 2d field (3d on history) at current lat on secondary histories:
!
      if (trim(name)=='ED1'.and.ied1 > 0) then
        call addfld(diags(ied1)%short_name,diags(ied1)%long_name,
     |    diags(ied1)%units,ed12,'mlon',mlon0,mlon1,diags(ied1)%levels,
     |    mlev0,mlev1,lat)
      elseif (trim(name)=='ED2'.and.ied2 > 0) then
        call addfld(diags(ied2)%short_name,diags(ied2)%long_name,
     |    diags(ied2)%units,ed12,'mlon',mlon0,mlon1,diags(ied2)%levels,
     |    mlev0,mlev1,lat)
      endif
      end subroutine mkdiag_ED12
!-----------------------------------------------------------------------
      subroutine mkdiag_HNMF2(name,ht,ne,lev0,lev1,lon0,lon1,lat)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: ht,ne
!
! Local:
      real,dimension(lon0:lon1) :: hmf2,nmf2,fof2
      integer,save :: ixhmf2=-1,ixnmf2=-1,ixfof2=-1

      if (trim(name)=='HMF2'.and.ixhmf2==-1) ixhmf2 = checkf(name)
      if (trim(name)=='NMF2'.and.ixnmf2==-1) ixnmf2 = checkf(name)
      if (trim(name)=='FOF2'.and.ixfof2==-1) ixfof2 = checkf(name)
      if (ixhmf2<=0.and.ixnmf2<=0.and.ixfof2<=0) return
!
! Calculate hmf2,nmf2,fof2 from ne at current latitude:
      call hnmf2(ht,ne,hmf2,nmf2,fof2,lev0,lev1,lon0,lon1,lat) 
!
! Save 2d field on secondary histories:
      if (trim(name)=='HMF2'.and.ixhmf2 > 0) then
        hmf2 = hmf2*1.e-5 ! cm to km
        call addfld(diags(ixhmf2)%short_name,diags(ixhmf2)%long_name,
     |    diags(ixhmf2)%units,hmf2,'lon',lon0,lon1,'lat',lat,lat,0)
      elseif (trim(name)=='NMF2'.and.ixnmf2 > 0) then
        call addfld(diags(ixnmf2)%short_name,diags(ixnmf2)%long_name,
     |    diags(ixnmf2)%units,nmf2,'lon',lon0,lon1,'lat',lat,lat,0)
      elseif (trim(name)=='FOF2'.and.ixfof2 > 0) then
        call addfld(diags(ixfof2)%short_name,diags(ixfof2)%long_name,
     |    diags(ixfof2)%units,fof2,'lon',lon0,lon1,'lat',lat,lat,0)
      endif
      end subroutine mkdiag_HNMF2
!-----------------------------------------------------------------------
      subroutine hnmf2(ht,ne,hmf2,nmf2,fof2,lev0,lev1,lon0,lon1,lat)
!
! Find the NmF2, foF2, and hmF2 based on the electron profile.
! btf: This is taken from an idl procedure provided by Liying.
! fof2 added April, 2012.
!
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: ht,ne
      real,intent(out),dimension(lon0:lon1) :: hmf2,nmf2,fof2
!
! Local:
      integer :: k,kx,i
      real :: h(3),p(3),h12,h22,h32,deltx,atx,ax,btx,bx,ctx,cx
!
      do i=lon0,lon1
        kx = 0
        kloop: do k = lev1-1,2,-1
          if (ne(k,i) >= ne(k-1,i).and.ne(k,i) >= ne(k+1,i)) then
            kx = k
            exit kloop
          endif
        enddo kloop
        if (kx==0) then
          write(6,"('>>> hnmf2: lat=',i4,' could not find kx -- ',
     |      'setting hmf2=nmf2=fof2=0')") lat
          hmf2(:) = 0. ; nmf2(:) = 0. ; fof2(:) = 0.
          return
        endif
        h = (/ht(kx-1,i),ht(kx,i),ht(kx+1,i)/)
        p = (/ne(kx-1,i),ne(kx,i),ne(kx+1,i)/)
        h12 = h(1)*h(1)
        h22 = h(2)*h(2)
        h32 = h(3)*h(3)
        deltx=h12*h(2)+h22*h(3)+h32*h(1)-h32*h(2)-h12*h(3)-h22*h(1)
        atx=p(1)*h(2)+p(2)*h(3)+p(3)*h(1)-h(2)*p(3)-h(3)*p(1)-h(1)*p(2)
!
! Periodic points are not set in ht and ne, so they are zero
! at i=1,2. Setting hmf2and nmf2 zero here is ok, since they
! are saved on histories at i=3,nlon. Without this check, code 
! will stop here w/ divide by zero deltx or zero atx when debug 
! flags are set.
!
        if (deltx == 0..or.atx == 0.) then
          hmf2(:) = 0. ; nmf2(:) = 0. ; fof2(:) = 0.
          cycle
        endif
        ax=atx/deltx
        btx=h12*p(2)+h22*p(3)+h32*p(1)-h32*p(2)-h12*p(3)-h22*p(1)
        bx=btx/deltx
        ctx=h12*h(2)*p(3)+h22*h(3)*p(1)+h32*h(1)*p(2)-h32*h(2)*p(1)-
     |    h12*h(3)*p(2)-h22*h(1)*p(3)
        cx=ctx/deltx
        hmf2(i)=-(bx/(2.*ax))
        nmf2(i)=-((bx*bx-4.*ax*cx)/(4.*ax))
        fof2(i)=sqrt(nmf2(i)*1.e6/1.24e10)
      enddo ! i=lon0,lon1
!     write(6,"('hnmf2: lat=',i4,' hmf2 min,max=',2e12.4,
!    |  ' nmf2 min,max=',2e12.4,' fof2 min,max=',2e12.4)") 
!    |  lat,minval(hmf2),maxval(hmf2),minval(nmf2),maxval(nmf2),
!    |  minval(fof2),maxval(fof2)
      end subroutine hnmf2
!-----------------------------------------------------------------------
!
! This is Liying's original IDL procedure hnmf2.pro. This was converted
! to fortran in subroutine hnmf2 above.
!
! PRO hnmf2, X, Y, XOUX, YOUX
! ;
! ; Find the NmF2 and hmF2 based on the electron profile
! ;
! ; input
! ;     X: height
! ;     Y: Ne profile
! ;
! ; output
! ;     XOUX: hmF2
! ;     YOUX: NmF2
! ;
! nz=n_elements(x)
! for k=nz-2,1,-1 do begin
!   if ((y[k] ge y[k-1]) and (y[k] ge y[k+1])) then break
! endfor
! h=[x[k-1],x[k],x[k+1]]
! p=[y[k-1],y[k],y[k+1]]
! 
! h12=h(0)*h(0)
! h22=h(1)*h(1)
! h32=h(2)*h(2)
! DELTX=h12*h(1)+h22*h(2)+h32*h(0)-h32*h(1)-h12*h(2)-h22*h(0)
! ATX=p(0)*h(1)+p(1)*h(2)+p(2)*h(0)-h(1)*p(2)-h(2)*p(0)-h(0)*p(1)
! AX=ATX/DELTX
! BTX=h12*p(1)+h22*p(2)+h32*p(0)-h32*p(1)-h12*p(2)-h22*p(0)
! BX=BTX/DELTX
! CTX=h12*h(1)*p(2)+h22*h(2)*p(0)+h32*h(0)*p(1)-h32*h(1)*p(0)-h12*h(2)*p(1)-h22*h(0)*p(2)
! CX=CTX/DELTX
! XOUX=-(BX/(2.*AX))
! YOUX=-((BX*BX-4.*AX*CX)/(4.*AX))
!
! RETURN
! END
!-----------------------------------------------------------------------
      subroutine mkdiag_SCHT(name,zcm,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate pressure scale height from geopotential Z (called from addiag.F)
! Note Z is input at interfaces, and pzps is output at midpoints
!   (as per conversation w/ Art)
!
!     diags(n)%short_name = 'SCHT'
!     diags(n)%long_name  = 'Pressure Scale Height'
!     diags(n)%units      = 'km'
!     diags(n)%levels     = 'lev'
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: zcm(lev0:lev1,lon0:lon1,lat0:lat1) ! geopotential (cm)
!
! Local:
      integer :: k,i,j
      real :: pzps(lev0:lev1,lon0:lon1)
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Take delta Z:
      do j=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            pzps(k,i) = zcm(k+1,i,j)-zcm(k,i,j)
          enddo
          pzps(lev1,i) = pzps(lev1-1,i)
!
! Generic for dlev 0.5 or 0.25 resolution:
          pzps(:,i) = pzps(:,i)/dlev
        enddo ! i=lon0,lon1
        pzps = pzps*1.e-5 ! cm to km
!
! Save on secondary history:
        call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |    diags(ix)%units,pzps,diags(ix)%levels,lev0,lev1,
     |    'lon',lon0,lon1,j)
      enddo ! j=lat0,lat1

      end subroutine mkdiag_SCHT
!-----------------------------------------------------------------------
      subroutine mkdiag_MU_M(name,mu_m,lev0,lev1,lon0,lon1,lat)
!
! Molecular viscosity coefficient (mu_m passed from cpktkm.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: mu_m
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,mu_m,diags(ix)%levels,lev0,lev1,
     |  'lon',lon0,lon1,lat)

      end subroutine mkdiag_MU_M
!-----------------------------------------------------------------------
      subroutine mkdiag_JQR(name,jqr,lon0,lon1,lat0,lat1)
!
! Upward current density at top of ionospheric shell 
! (jqr passed from current.F90 with current_kq==1 in dynamo.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lat0,lat1,lon0,lon1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: jqr
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,jqr,'mlon',lon0,lon1,'mlat',lat0,lat1,0)

      end subroutine mkdiag_JQR
!-----------------------------------------------------------------------
      subroutine mkdiag_KQLAM(name,kqlam,lon0,lon1,lat0,lat1)
!
! magn. northward height-integrated current density
! (kqlam passed from current.F90 with current_kq==1 in dynamo.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lat0,lat1,lon0,lon1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: kqlam
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,kqlam,'mlon',lon0,lon1,'mlat',lat0,lat1,0)

      end subroutine mkdiag_KQLAM
!-----------------------------------------------------------------------
      subroutine mkdiag_KQPHI(name,kqphi,lon0,lon1,lat0,lat1)
!
! magn. eastward height-integrated current density 
! (kqphi passed from current.F90 with current_kq==1 in dynamo.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lat0,lat1,lon0,lon1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: kqphi
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,kqphi,'mlon',lon0,lon1,'mlat',lat0,lat1,0)

      end subroutine mkdiag_KQPHI
!-----------------------------------------------------------------------
      subroutine mkdiag_JE13D(name,je13d,lon0,lon1,lev0,lev1,lat)
!
! Eastward current density (geomagnetic grid)
! (je13d passed from current.F90 with current_kq==1 in dynamo.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat
      real,intent(in),dimension(lon0:lon1,lev0:lev1) :: je13d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,je13d,'mlon',lon0,lon1,'mlev',lev0,lev1,lat)

      end subroutine mkdiag_JE13D
!-----------------------------------------------------------------------
      subroutine mkdiag_JE23D(name,je23d,lon0,lon1,lev0,lev1,lat)
!
! Equator-/downward current density (geomagnetic grid)
! (je23d passed from current.F90 with current_kq==1 in dynamo.F):
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat
      real,intent(in),dimension(lon0:lon1,lev0:lev1) :: je23d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,je23d,'mlon',lon0,lon1,'mlev',lev0,lev1,lat)

      end subroutine mkdiag_JE23D
!-----------------------------------------------------------------------
      subroutine mkdiag_WN(name,omega,zcm,lev0,lev1,lon0,lon1,lat)
!
! Neutral Vertical Wind on interfaces, from vertical motion OMEGA and 
! scale height. Scale height pzps is calculated from input geopotential z (cm).
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: omega,zcm
!
! Local:
      integer :: i,k
      real,dimension(lev0:lev1,lon0:lon1) :: wn
      real,dimension(lev0:lev1) :: pzps,omega1
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! Calculate scale height pzps:
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          pzps(k) = (zcm(k+1,i)-zcm(k-1,i))/(2.*dlev)
        enddo
        pzps(lev0) = (zcm(lev0+1,i)-zcm(lev0,i))/dlev
        pzps(lev1) = pzps(lev1-1)
!
        omega1(:) = omega(:,i)
        omega1(lev1) = omega1(lev1-1)
!
! Output vertical wind (cm):
        wn(:,i) = omega1(:)*pzps(:)
      enddo ! i=lon0,lon1

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,wn,diags(ix)%levels,lev0,lev1,
     |  'lon',lon0,lon1,lat)

      end subroutine mkdiag_WN
!-----------------------------------------------------------------------
      subroutine mkdiag_O_N2(name,o1,o2,n2,lev0,lev1,lon0,lon1,lat)
!
! Calculate O/N2 ratio from o2 and o (mmr). 
! N2 is passed in (3d array in fields.F), where n2=1-o2-o-he
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: o1,o2,n2
!
! Local:
      real,dimension(lev0:lev1,lon0:lon1) :: o_n2
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! O/N2 ratio:
      o_n2 = o1/n2

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,o_n2,'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine mkdiag_O_N2
!-----------------------------------------------------------------------
      subroutine mkdiag_N2(name,n2,lev0,lev1,lon0,lon1,lat)
!
! Save updated N2 (3d field in fields.F) from comp.F:
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: n2
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,n2,'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine mkdiag_N2
!-----------------------------------------------------------------------
      subroutine mkdiag_QJOULE_INTEG(name,qji_tn,lev0,lev1,lon0,lon1,
     |  lat)
      use cons_module,only: p0,grav
      use params_module,only: zpint
!
! Calculate height-integrated Joule heating (called from qjoule.F)
! This method is adapted from ncl code provided by Astrid (7/20/11)
!
! Args
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: qji_tn
!
! Local:
      integer :: k,i
      real,dimension(lon0:lon1) :: qji_integ
      real,dimension(lev0:lev1,lon0:lon1) :: qj
      real :: myp0,mygrav
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return
!
! First integrate to get MKS units W/m^2:
! (If you want these units, comment out the below conversion to CGS)
!
      mygrav = grav*.01     ! cm/s^2 to m/s^2
      myp0 = p0*1.e-3*100.  ! to Pa
      qj = qji_tn*.0001     ! ergs/g/s to W/kg 10^(-7)*10^3
      
      qji_integ = 0.
      do i=lon0,lon1
        do k=lev0,lev1-1
          qji_integ(i) = qji_integ(i) + myp0/mygrav*exp(-zpint(k))*
     |      qj(k,i)*dlev
        enddo
      enddo
!
! Output in CGS units, to be consistent w/ the model:
! (note that 1 erg/cm^2/s == 1 mW/m^2)
      qji_integ = qji_integ*1000. ! W/m^2 to erg/cm^2/s
!
! Save 2d field on secondary history:
      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,qji_integ,'lon',lon0,lon1,'lat',lat,lat,0)

      end subroutine mkdiag_QJOULE_INTEG
!-----------------------------------------------------------------------
      subroutine mkdiag_PHIM2D(name,phim2d,mlon0,mlon1,mlat0,mlat1)
!
! 2d Electric Potential on magnetic grid (called from dynamo with
! global mag dimensions)
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: mlon0,mlon1,mlat0,mlat1
      real,intent(in),dimension(mlon0:mlon1,mlat0:mlat1) :: phim2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,phim2d,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

      end subroutine mkdiag_PHIM2D
!-----------------------------------------------------------------------
      subroutine mkdiag_ZGMID(name,zgint,lev0,lev1,lon0,lon1,lat)
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,intent(in),dimension(lev0:lev1,lon0:lon1) :: 
     |  zgint ! geometric height at interfaces
!
! Local:
      integer,save :: ix=0,firstcall=0,k
      real,dimension(lev0:lev1,lon0:lon1) :: zgmid

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      do k=lev0,lev1-1
        zgmid(k,:) = 0.5*(zgint(k,:)+zgint(k+1,:))
      enddo
      zgmid(lev1,:) = 2.0*zgmid(lev1-1,:)-zgmid(lev1-2,:)

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,zgmid,'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine mkdiag_ZGMID
!-----------------------------------------------------------------------
      subroutine mkdiag_CUSP(name,cusp2d,lon0,lon1,lat0,lat1)
      use aurora_module,only: ec
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: cusp2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,cusp2d*ec,'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine mkdiag_CUSP
!-----------------------------------------------------------------------
      subroutine mkdiag_DRIZZLE(name,drzl2d,lon0,lon1,lat0,lat1)
      use aurora_module,only: ed
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: drzl2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,drzl2d*ed,'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine mkdiag_DRIZZLE
!-----------------------------------------------------------------------
      subroutine mkdiag_ALFA(name,alfa2d,lon0,lon1,lat0,lat1)
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: alfa2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,alfa2d,'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine mkdiag_ALFA
!-----------------------------------------------------------------------
      subroutine mkdiag_NFLUX(name,nflx2d,lon0,lon1,lat0,lat1)
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: nflx2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,nflx2d,'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine mkdiag_NFLUX
!-----------------------------------------------------------------------
      subroutine mkdiag_EFLUX(name,eflux2d,lon0,lon1,lat0,lat1)
!
! Interpolate geometric height on interfaces (input zgint) to midpoints.
!
! Args:
      character(len=*),intent(in) :: name
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in),dimension(lon0:lon1,lat0:lat1) :: eflux2d
!
! Local:
      integer,save :: ix=0,firstcall=0

      firstcall = firstcall+1
      if (firstcall == 1) ix = checkf(name)
      if (ix==0) return

      call addfld(diags(ix)%short_name,diags(ix)%long_name,
     |  diags(ix)%units,eflux2d,'lon',lon0,lon1,'lat',lat0,lat1,0)

      end subroutine mkdiag_EFLUX
!-----------------------------------------------------------------------
!
! End diags module:
      end module diags_module
