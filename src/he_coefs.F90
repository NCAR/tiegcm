module he_coefs_module
  implicit none

  integer :: nlat_he,nlon_he,ntrigs
  integer,dimension(13) :: ifax
  real,dimension(:),allocatable :: trigs
  real,dimension(:),allocatable :: glat_he,glon_he
  real,dimension(:,:,:),allocatable :: pmn,zmn

  contains
!-----------------------------------------------------------------------
  subroutine init_he_coefs

    use input_module,only: he_coefs_ncfile
    use netcdf,only: nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_strerror,nf90_nowrite,nf90_noerr

    integer :: stat,ncid,dimid,varid,length,i
    real :: res
    external :: shutdown

    stat = nf90_open(trim(he_coefs_ncfile),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at opening '//trim(he_coefs_ncfile))

    stat = nf90_inq_dimid(ncid,'lat1',dimid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat1 ID')

    stat = nf90_inquire_dimension(ncid,dimid,len=nlat_he)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat1 length')

    stat = nf90_inq_dimid(ncid,'lat2',dimid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat2 ID')

    stat = nf90_inquire_dimension(ncid,dimid,len=length)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat2 length')

    if (length /= nlat_he-1) call shutdown('Invalid Helium coefficient file')

    stat = nf90_inq_dimid(ncid,'lat3',dimid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat3 ID')

    stat = nf90_inquire_dimension(ncid,dimid,len=length)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat3 length')

    if (length /= nlat_he) call shutdown('Invalid Helium coefficient file')

    allocate(pmn(nlat_he,0:nlat_he-2,0:nlat_he-1))
    allocate(zmn(nlat_he,0:nlat_he-2,0:nlat_he-1))

    stat = nf90_inq_varid(ncid,'pmn',varid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring variable pmn ID')

    stat = nf90_get_var(ncid,varid,pmn)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at getting variable pmn')

    stat = nf90_inq_varid(ncid,'zmn',varid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at inquiring variable zmn ID')

    stat = nf90_get_var(ncid,varid,zmn)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at getting variable zmn')

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) &
      call shutdown(trim(nf90_strerror(stat))//' at closing '//trim(he_coefs_ncfile))

    nlon_he = nlat_he*2
    res = 180./nlat_he

    allocate(glat_he(nlat_he))
    allocate(glon_he(nlon_he))

    glat_he = (/(-90.-res/2.+i*res,i=1,nlat_he)/)
    glon_he = (/(-180.+(i-1)*res,i=1,nlon_he)/)

    ntrigs = 3*nlon_he/2+1
    allocate(trigs(ntrigs))

  end subroutine init_he_coefs
!-----------------------------------------------------------------------
end module he_coefs_module
