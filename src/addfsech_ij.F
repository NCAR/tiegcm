!
      subroutine addfsech_ij(name,long_name,units,f2d,lon0,lon1,
     |  lat0,lat1)
!
! Add a 2d (lon,lat) field (geo, mag, or magphr) to secondary histories. 
! This will define fsech2d(ix)%data(lon0:lon1,lat0:lat1) (or fsechmag2d)
!   from the input field f2d.
! The lon and lat ranges may be geographic, magnetic, or magphr.
!
      use fields_module,only: fsech2d,fsechmag2d,fsechmagphr2d,
     |  f4d,nf4d,shortname_len
      use hist_module,only: isechist,modeltime
      use init_module,only: istep
      use params_module,only: mxfsech
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,lat0,lat1
      real,intent(in) :: f2d(lon0:lon1,lat0:lat1)
!
! Local:
      integer :: i,ix,k,nlev_used,iseries_p,iseries_s,ixgeo2d,ixmag2d,
     |  ixmagphr
      real :: fmin,fmax
      logical :: time2write,wrprim,wrsech,save_prim,save_sech,
     |  newseries_prim,newseries_sech
!
! External:
      integer,external :: strloc
      logical,external :: wrhist
!
! Return silently if not writing secondary histories:
      if (isechist==0) return
!
! Save only if we are writing a secondary history this step:
      time2write = wrhist(istep,modeltime,
     |  wrprim, save_prim, newseries_prim, iseries_p,
     |  wrsech, save_sech, newseries_sech, iseries_s)
      if (.not.wrsech) return
!
!     write(6,"('addfsech_ij: lon0,1=',2i3,' lat0,1=',2i3,' istep=',i6,
!    |  ' modeltime=',4i4,' name=',a)") lon0,lon1,lat0,lat1,istep,
!    |  modeltime,name
!
! Check if field is a prognostic (i.e., a 4-d field):
      ix = strloc(f4d%short_name,nf4d,name)
      if (ix > 0) then
        write(6,"('>>> addfsech_ij: field ',a,' is a prognostic',
     |    ' -- returning.')") trim(name)
        return
      endif
!
! Check for geographic, magnetic, or magphr:
      ixgeo2d  = strloc(fsech2d      %short_name,mxfsech,name)
      ixmag2d  = strloc(fsechmag2d   %short_name,mxfsech,name)
      ixmagphr = strloc(fsechmagphr2d%short_name,mxfsech,name)
!
! Call appropriate subroutine:
      if (ixgeo2d > 0) then ! is geographic2d
        call addfsechij_geo2d(name,long_name,units,f2d,ixgeo2d,lon0,
     |    lon1,lat0,lat1)
      elseif (ixmag2d > 0) then ! is geomagnetic2d
        call addfsechij_mag2d(name,long_name,units,f2d,ixmag2d,lon0,
     |    lon1,lat0,lat1)
      elseif (ixmagphr > 0) then ! is magnetosphere2d
        call addfsechij_magphr(name,long_name,units,f2d,ixmagphr,lon0,
     |    lon1,lat0,lat1)
      else
!       write(6,"('>>> addfsech_ij: field ',a,
!    |    ' not found in SECFGEO2D, SECFMAG2D, or SECFMAGPHR.')")
!    |    trim(name)
        return
      endif
      end subroutine addfsech_ij
!-----------------------------------------------------------------------
      subroutine addfsechij_geo2d(name,long_name,units,f2d,ix,lon0,lon1,
     |  lat0,lat1)
!
! Add geographic f2d(lon,lat) to secondary histories, i.e., 
!   define fsech2d(ix)%data(lon0:lon1,lat0:lat1).
!
      use fields_module,only: fsech2d
      use params_module,only: nlonp4,nlat,nlon
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,lat0,lat1,ix
      real,intent(in) :: f2d(lon0:lon1,lat0:lat1)
!
! Check dimensions:
      if (lon0 < 1 .or. lon1 > nlonp4 .or. lon0 > lon1) then
        write(6,"('>>> addfsechij_geo2d: lon0 must be >= 1, lon1 ',
     |    'must be <= nlonp4, and lon0 must be <= lon1')")
        write(6,"('    lon0=',i3,' lon1=',i3,' nlonp4=',i3)")
     |    lon0,lon1,nlonp4
        call shutdown('addfsechij_geo2d')
      endif
      if (lat0 < 1 .or. lat1 > nlat .or. lat0 > lat1) then
        write(6,"('>>> addfsechij_geo2d: lat0 must be >= 1, lat1 ',
     |    'must be <= nlat, and lat0 must be <= lat1')")
        write(6,"('    lat0=',i3,' lat1=',i3,' nlat=',i3)")
     |    lat0,lat1,nlat
        call shutdown('addfsechij_geo2d')
      endif
!
! fsech2d was allocated in set_fsech2d (flds_mod.f):
! Double check that this pointer is ok:
!
      if (.not.associated(fsech2d(ix)%data)) then
        write(6,"('>>> addfsechij_geo2d: fsech2d(ix)%data data not ',
     |    'associated. ix=',i3)") ix
        call shutdown('addfsechij_geo2d')
      endif
 !
! Define fsech(ix)%data for current latitude:
! Save 2-d field:
      fsech2d(ix)%data(lon0:lon1,lat0:lat1) = f2d(lon0:lon1,lat0:lat1)
!
      fsech2d(ix)%long_name = long_name
      fsech2d(ix)%units     = units
      end subroutine addfsechij_geo2d
!-----------------------------------------------------------------------
      subroutine addfsechij_mag2d(name,long_name,units,f2d,ix,lon0,
     |  lon1,lat0,lat1)
!
! Add geomagnetic f2d(lon,lat) to secondary histories, i.e., 
!   define fsechmag2d(ix)%data(lon0:lon1,lat0:lat1).
!
      use fields_module,only: fsechmag2d
      use params_module,only: nmlonp1,nmlat,nmlon
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,ix,lat0,lat1
      real,intent(in) :: f2d(lon0:lon1,lat0:lat1)
!
! Check dimensions:
      if (lon0 < 1 .or. lon1 > nmlonp1 .or. lon0 > lon1) then
        write(6,"('>>> addfsechij_mag2d: lon0 must be >= 1, lon1 must ',
     |    'be ',' <= nmlonp1, and lon0 must be <= lon1')")
        write(6,"('    lon0=',i3,' lon1=',i3,' nmlonp1=',i3)")
     |    lon0,lon1,nmlonp1
        call shutdown('addfsechij_mag2d')
      endif
      if (lat0 < 1 .or. lat1 > nmlat .or. lat0 > lat1) then
        write(6,"('>>> addfsechij_mag2d: lat0 must be >= 1, lat1 must ',
     |    'be ',' <= nmlat, and lat0 must be <= lat1')")
        write(6,"('    lat0=',i3,' lat1=',i3,' nmlat=',i3)")
     |    lat0,lat1,nmlat
        call shutdown('addfsechij_mag2d')
      endif
!
! fsechmag2d was allocated in set_fsechmag2d (flds_mod.f):
!   allocate(fsechmag2d(nfsech)%data(nmlonp1,nmlat),
! Double check that this pointer is ok:
!
      if (.not.associated(fsechmag2d(ix)%data)) then
        write(6,"('>>> addfsechij_mag2d: fsechmag2d(ix)%data data not ',
     |    'associated. ix=',i3)") ix
        call shutdown('addfsechij_mag2d')
      endif
!
! Define fsechmag2d(ix)%data for current latitude:
!
      fsechmag2d(ix)%data(lon0:lon1,lat0:lat1)=f2d(lon0:lon1,lat0:lat1)
!
      fsechmag2d(ix)%long_name = long_name
      fsechmag2d(ix)%units     = units
      end subroutine addfsechij_mag2d
!-----------------------------------------------------------------------
      subroutine addfsechij_magphr(name,long_name,units,f2d,ix,lon0,
     |  lon1,lat0,lat1)
!
! Add magnetospheric f2d(lon,lat) to secondary histories, i.e., 
!   define fsechmagphr2d(ix)%data(lon0:lon1,lat0:lat1).
!
      use fields_module,only: fsechmagphr2d
      use params_module,only: nmagphrlon,nmagphrlat
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,ix,lat0,lat1
      real,intent(in) :: f2d(lon0:lon1,lat0:lat1)
!
! Check dimensions:
      if (lon0 < 1 .or. lon1 > nmagphrlon .or. lon0 > lon1) then
        write(6,"('>>> addfsechij_magphr: lon0 must be >= 1, lon1 ',
     |    'must be ',' <= nmagphrlon, and lon0 must be <= lon1')")
        write(6,"('    lon0=',i3,' lon1=',i3,' nmagphrlon=',i3)")
     |    lon0,lon1,nmagphrlon
        call shutdown('addfsechij_magphr')
      endif
      if (lat0 < 1 .or. lat1 > nmagphrlat .or. lat0 > lat1) then
        write(6,"('>>> addfsechij_magphr: lat0 must be >= 1, lat1 ',
     |    'must be ',' <= nmagphrlat, and lat0 must be <= lat1')")
        write(6,"('    lat0=',i3,' lat1=',i3,' nmagphrlat=',i3)")
     |    lat0,lat1,nmagphrlat
        call shutdown('addfsechij_magphr')
      endif
!
! fsechmagphr2d was allocated in set_fsechmagphr2d (flds_mod.f):
!   allocate(fsechmagphr2d(nfsech)%data(nmagphrlon,nmagphrlat),
! Double check that this pointer is ok:
!
      if (.not.associated(fsechmagphr2d(ix)%data)) then
        write(6,"('>>> addfsechij_magphr: fsechmagphr(ix)%data data',
     |    ' not associated. ix=',i3)") ix
        call shutdown('addfsechij_magphr')
      endif
!
! Define fsechmagphr2d(ix)%data for current latitude:
!
      fsechmagphr2d(ix)%data(lon0:lon1,lat0:lat1) =
     |   f2d(lon0:lon1,lat0:lat1)
!
      fsechmagphr2d(ix)%long_name = long_name
      fsechmagphr2d(ix)%units     = units
      end subroutine addfsechij_magphr

