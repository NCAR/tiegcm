!
      subroutine addfsech(name,long_name,units,f2d,lon0,lon1,levdim,
     |  nlevreq,lat)
!
! Add a field to secondary histories.
! This routine is user-callable from inside the main latitude loop.
! It adds a 2d latitude slice of a diagnostic field to fsech(i)%data,
! which will be written to the secondary histories.
!
      use fields_module,only: fsech,fsechmag,f4d,nf4d,shortname_len
      use hist_module,only: isechist,modeltime
      use init_module,only: istep
      use params_module,only: mxfsech
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,levdim,nlevreq,lat
      real,intent(in) :: f2d(levdim,lon0:lon1)
!
! Local:
      integer :: i,ix,ixmag,k,nlev_used,istat,iseries_p,iseries_s
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
!     write(6,"('addfsech: lat=',i2,' lon0,1=',2i3,' istep=',i6,
!    |  ' modeltime=',4i4,' name=',a)") lat,lon0,lon1,istep,modeltime,
!    |  name
!
! Check if field is a prognostic (i.e., a 4-d field):
      ix = strloc(f4d%short_name,nf4d,name) 
      if (ix > 0) then
        write(6,"('>>> addfsech: field ',a,' is a prognostic',
     |    ' -- returning.')") trim(name)
        return
      endif
!
! Check that field was requested on secondary histories:
! (i.e., was on either SECFLDS or SECFMAG inputs)
! Check geographic fields:
      ixmag = 0
      ix = strloc(fsech%short_name,mxfsech,name) 
      if (ix == 0) then ! was not in SECFLDS (not geographic)
!
! Check magnetic fields:
        ixmag = strloc(fsechmag%short_name,mxfsech,name) 
        if (ixmag==0) then
!         if (lat==2) write(6,"('>>> addfsech: field ',a,
!    |      ' not found in requested secondary history field ',
!    |      'names.')") trim(name)
          return
        endif
      endif
!
! Call for geographic or magnetic:
      if (ix > 0) then
!       write(6,"('Calling addfsech_geo for field ',a,' at lat=',i3,
!    |    ' ix=',i3)") name,lat,ix
        call addfsech_geo(name,long_name,units,f2d,ix,lon0,lon1,levdim,
     |    nlevreq,lat)
      else
!       write(6,"('Calling addfsech_mag for field ',a,' at lat=',i3)")
!    |    name,lat
        call addfsech_mag(name,long_name,units,f2d,ixmag,lon0,lon1,
     |    levdim,nlevreq,lat)
      endif
      return
      end subroutine addfsech
!-----------------------------------------------------------------------
      subroutine addfsech_geo(name,long_name,units,f2d,ix,lon0,lon1,
     |  levdim,nlevreq,lat)
!
! Add a field on the geographic grid to secondary histories.
! This adds a 2d latitude slice of a diagnostic field to fsech(i)%data,
! which will be written to the secondary histories.
!
      use fields_module,only: fsech
      use params_module,only: nlonp4,nlevp1,nlat,spval
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,levdim,nlevreq,lat,ix
      real,intent(in) :: f2d(levdim,lon0:lon1)
!
! Local:
      integer :: k,nlev_used
      real :: fmin,fmax
!
! Check dimensions:
      if (lon0 < 1 .or. lon1 > nlonp4 .or. lon0 > lon1) then
        write(6,"('>>> addfsech_geo: lon0 must be >= 1, lon1 must be ',
     |    ' <= nlonp4, and lon0 must be <= lon1')")
        write(6,"('    lon0=',i3,' lon1=',i3,' nlonp4=',i3)")
     |    lon0,lon1,nlonp4
        call shutdown('addfsech_geo')
      endif
      if (levdim /= nlevp1) then
        write(6,"('>>> addfsech_geo: levdim must equal nlevp1: levdim=',
     |    i4,' nlevp1=',i4)") levdim,nlevp1
        call shutdown('addfsech_geo')
      endif
      if (lat <= 0 .or. lat > nlat) then
        write(6,"('>>> addfsech_geo: bad lat=',i4,' (nlat=',i3,')')")
     |    lat,nlat
        call shutdown('addfsech_geo')
      endif
!
! nlev must be <= levdim:
      nlev_used = nlevreq
      if (nlevreq > levdim) then 
        write(6,"('>>> addfsech_geo: requested nlev must be <= levdim:',
     |    ' nlevreq=',i4,' levdim=',i4,' Will use nlevreq==levdim')") 
     |  nlev_used,levdim
        nlev_used = levdim
      endif
!
! fsech was allocated in set_fsech (flds_mod.f):
! Double check that this pointer is ok:
!
      if (.not.associated(fsech(ix)%data)) then
        write(6,"('>>> addfsech_geo: fsech(ix)%data data not ',
     |    'associated. ix=',i3)") ix 
        call shutdown('addfsech_geo')
      endif
!
! Define fsech(ix)%data for current latitude:
! (if nlev < levdim, pad top levels nlev+1,levdim with spval.
!
      if (nlev_used == levdim) then
        fsech(ix)%data(:,lon0:lon1,lat) = f2d(:,lon0:lon1) 
      else ! nlev < levdim
        do k=1,nlev_used
          fsech(ix)%data(k,lon0:lon1,lat) = f2d(k,lon0:lon1)
        enddo
        do k=nlev_used+1,levdim
          fsech(ix)%data(k,lon0:lon1,lat) = spval
        enddo
      endif

!     call fminmax(fsech(ix)%data(1:nlev_used,lon0:lon1,lat),
!    |  nlev_used*(lon1-lon0+1),fmin,fmax)
!     write(6,"('addfsech_geo: lat=',i3,' lon0,1=',2i3,' field ',a,
!    |  ' min,max=',2e12.4)") lat,lon0,lon1,fsech(ix)%short_name(1:4),
!    |  fmin,fmax
!
      fsech(ix)%long_name = long_name
      fsech(ix)%units = units
!
! Set task0_only: If lon0-lon1 covers global longitude domain, then
! assume the field is being defined by root task only. If task0_only
! is true, then gather2root will not be done on this field. (task0_only 
! was initialized to false in sub set_fsech (fields.F)).
!
      if (lon0==1.and.lon1 >= nlonp4-3) fsech(ix)%task0_only = .true.
      end subroutine addfsech_geo
!-----------------------------------------------------------------------
      subroutine addfsech_mag(name,long_name,units,f2d,ix,lon0,lon1,
     |  levdim,nlevreq,lat)
!
! Add a field on the magnetic grid field to secondary histories.
! This adds a 2d latitude slice of a diagnostic field to fsechmag(i)%data,
! which will be written to the secondary histories.
!
      use fields_module,only: fsechmag
      use params_module,only: nmlonp1,nmlat,nlevp1,spval
      implicit none
!
! Args:
      character(len=*),intent(in) :: name,long_name,units
      integer,intent(in) :: lon0,lon1,levdim,nlevreq,lat,ix
      real,intent(in) :: f2d(levdim,lon0:lon1)
!
! Local:
      integer :: k,nlev_used
      real :: fmin,fmax
!
! Check dimensions:
      if (lon0 < 1 .or. lon1 > nmlonp1 .or. lon0 > lon1) then
        write(6,"('>>> addfsech_mag: lon0 must be >= 1, lon1 must be ',
     |    ' <= nmlonp1, and lon0 must be <= lon1')")
        write(6,"('    lon0=',i3,' lon1=',i3,' nmlonp1=',i3)")
     |    lon0,lon1,nmlonp1
        call shutdown('addfsech_mag')
      endif
      if (levdim /= nlevp1+3) then
        write(6,"('>>> addfsech_mag: levdim must equal nlevp1+3: ',
     |    'levdim=',i4,' nlevp1+3=',i4)") levdim,nlevp1+3
        call shutdown('addfsech_mag')
      endif
      if (lat <= 0 .or. lat > nmlat) then
        write(6,"('>>> addfsech_mag: bad lat=',i4,' (nmlat=',i3,')')")
     |    lat,nmlat
        call shutdown('addfsech_mag')
      endif
!
! nlev must be <= levdim:
      nlev_used = nlevreq
      if (nlevreq > levdim) then 
        write(6,"('>>> addfsech_mag: requested nlev must be <= levdim:',
     |    ' nlevreq=',i4,' levdim=',i4,' Will use nlev==levdim')") 
     |    nlevreq,levdim
        nlev_used = levdim
      endif
!
! fsechmag was allocated in set_fsechmag (flds_mod.f):
!   allocate(fsechmag(nfsech)%data(nmlonp1,nmlat,nlevp1+3),
! Double check that this pointer is ok:
!
      if (.not.associated(fsechmag(ix)%data)) then
        write(6,"('>>> addfsech_mag: fsechmag(ix)%data data not ',
     |    'associated. ix=',i3)") ix 
        call shutdown('addfsech_mag')
      endif
!
! Define fsechmag(ix)%data for current latitude:
! (if nlev < levdim, pad top levels (nlev+1 to levdim) with spval.
!
      do k=1,nlev_used
        fsechmag(ix)%data(k,lon0:lon1,lat) = f2d(k,lon0:lon1) 
      enddo
      if (nlev_used < levdim) then 
        do k=nlev_used+1,levdim
          fsechmag(ix)%data(k,lon0:lon1,lat) = spval
        enddo
      endif
!
      fsechmag(ix)%long_name = long_name
      fsechmag(ix)%units = units
!
! Set task0_only: If lon0-lon1 covers global longitude domain, then
! assume the field is being defined by root task only. If task0_only
! is true, then gather2root will not be done on this field. (task0_only 
! was initialized to false in sub set_fsech (fields.F)).
!
      if (lon0==1.and.lon1 >= nmlonp1-1) 
     |  fsechmag(ix)%task0_only = .true.
      end subroutine addfsech_mag
