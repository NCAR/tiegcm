module nudge_module
! force the model fields near lower boundary with external fields

  implicit none

! if external fields contain lb, then the model lower boundary will use external fields
  character(len=*), dimension(*), parameter :: lb = (/'TN', 'UN', 'VN', 'Z '/)
  integer, parameter :: nlb = size(lb)

! whether the external field has a full longitude cycle
  logical :: wrap

  integer :: nlon, nlat, nlev, nilev, ntime, nfile, nf4d, maxlev, &
    ifile, itime, ncid, latbeg, latend, lonbeg, lonend, offbeg, offend, lonbeg1, lonend1
  logical, dimension(:), allocatable :: vcoord_is_midpoint
  integer, dimension(:), allocatable :: time, nt, t0, t1, varid
  real, dimension(:), allocatable :: lon, lat, lev, ilev

! whether this lower boundary field is to be used 
  logical, dimension(nlb) :: lbc_flag

! index of external fields in the model fields
  integer, dimension(:), allocatable :: f4d_idx

! vertical and horizontal relaxation factors
  real, dimension(:,:), allocatable :: vert_weight, hori_weight

  real, dimension(:,:,:,:), allocatable :: lbc
  real, dimension(:,:,:,:,:), allocatable :: f4d

  contains
!-----------------------------------------------------------------------
  subroutine check
! check the compliance of the external field

    use params_module, only: mxhvols, nlevp1, zibot, zpmid, zpint
    use input_module, only: mxlen_filename, nudge_ncpre, nudge_ncfile, nudge_ncpost, &
      nudge_flds, nudge_lbc, nudge_sponge, nudge_use_refdate, nudge_refdate, start_year, start_day
    use fields_module, only: f4d_int=>f4d
    use char_module, only: find_index
    use netcdf, only: nf90_open, nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
      nf90_get_var, nf90_inquire_variable, nf90_close, nf90_strerror, nf90_nowrite, nf90_noerr

    integer, parameter :: maxnt = 1440
    integer :: stat, dimid_lon, dimid_lat, dimid_lev, dimid_ilev, dimid_time, &
      varid_lon, varid_lat, varid_lev, varid_ilev, varid_time, varid_date, varid_datesec, &
      ifld, lb_idx, ndims, ik, it, yr, mn, dy, dt0, start_datenum
    real :: dlon
    character(len=mxlen_filename) :: filename
    integer, dimension(4) :: dimids
    integer, dimension(mxhvols, maxnt) :: t, dt, dts
    logical, external :: isclose
    integer, external :: to_doy, to_datenum
    external :: shutdown

    filename = trim(nudge_ncpre)//trim(nudge_ncfile(1))//trim(nudge_ncpost)

    stat = nf90_open(trim(filename), nf90_nowrite, ncid)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at opening '//trim(filename))

    stat = nf90_inq_dimid(ncid, 'lon', dimid_lon)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lon ID')

    stat = nf90_inquire_dimension(ncid, dimid_lon, len=nlon)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lon length')

    stat = nf90_inq_dimid(ncid, 'lat', dimid_lat)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat ID')

    stat = nf90_inquire_dimension(ncid, dimid_lat, len=nlat)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lat length')

    stat = nf90_inq_dimid(ncid, 'lev', dimid_lev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lev ID')

    stat = nf90_inquire_dimension(ncid, dimid_lev, len=nlev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension lev length')

    stat = nf90_inq_dimid(ncid, 'ilev', dimid_ilev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension ilev ID')

    stat = nf90_inquire_dimension(ncid, dimid_ilev, len=nilev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension ilev length')

    stat = nf90_inq_dimid(ncid, 'time', dimid_time)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension time ID')

    allocate(lon(nlon))
    allocate(lat(nlat))
    allocate(lev(nlev))
    allocate(ilev(nilev))

    stat = nf90_inq_varid(ncid, 'lon', varid_lon)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable lon ID')

    stat = nf90_get_var(ncid, varid_lon, lon)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable lon')

! if the external field covers the full longitude cycle,
! then model fields at lon=1,2 and lon=nlonp4-1,nlonp4 will be obtained from the external field
    dlon = (lon(nlon)-lon(1)) / (nlon-1)
    if (isclose(lon(1)+360-lon(nlon), dlon)) then
      wrap = .true.
    else
      wrap = .false.
    endif

    stat = nf90_inq_varid(ncid, 'lat', varid_lat)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable lat ID')

    stat = nf90_get_var(ncid, varid_lat, lat)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable lat')

    stat = nf90_inq_varid(ncid, 'lev', varid_lev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable lev ID')

    stat = nf90_get_var(ncid, varid_lev, lev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable lev')

    if (lev(1)>zibot .or. lev(nlev)<zibot) call shutdown('nudge_ncfile dimension lev must include model lbc')

    stat = nf90_inq_varid(ncid, 'ilev', varid_ilev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable ilev ID')

    stat = nf90_get_var(ncid, varid_ilev, ilev)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable ilev')

    if (ilev(1)>zibot .or. ilev(nilev)<zibot) call shutdown('nudge_ncfile dimension ilev must include model lbc')

    do ik = 1, nlevp1
      if (max(zpmid(ik), zpint(ik)) > zibot+nudge_sponge(2)) exit
    enddo
    maxlev = ik - 1

    if (nudge_sponge(1)*2>=max(lon(nlon)-lon(1), lat(nlat)-lat(1)) .or. &
      nudge_sponge(2)>max(lev(nlev)-lev(1), ilev(nilev)-ilev(1))) &
      call shutdown('nudge_sponge cannot exceed nudge_ncfile dimension range')

    nf4d = count(len_trim(nudge_flds) > 0)
    allocate(vcoord_is_midpoint(nf4d))
    allocate(varid(nf4d))
    allocate(f4d_idx(nf4d))
    lbc_flag = .false.
    do ifld = 1, nf4d
      f4d_idx(ifld) = find_index(nudge_flds(ifld), f4d_int%short_name)
      if (f4d_idx(ifld) == 0) call shutdown(trim(nudge_flds(ifld))//' is not a valid model field')

      lb_idx = find_index(nudge_flds(ifld), lb)
      if (nudge_lbc(ifld)) then
        if (lb_idx == 0) then
          call shutdown('nudge_flds '//nudge_flds(ifld)//' is not a lower boundary field')
        else
          lbc_flag(lb_idx) = .true.
        endif
      endif
    enddo

    do ifld = 1, nf4d
      stat = nf90_inq_varid(ncid, trim(nudge_flds(ifld)), varid(ifld))
      if (stat /= nf90_noerr) &
        call shutdown(trim(nf90_strerror(stat))//' at inquiring variable '//trim(nudge_flds(ifld))//' ID')

      stat = nf90_inquire_variable(ncid, varid(ifld), ndims=ndims, dimids=dimids)
      if (stat /= nf90_noerr) &
        call shutdown(trim(nf90_strerror(stat))//' at inquiring variable '//trim(nudge_flds(ifld))//' dimensions')

      if (ndims /= 4) call shutdown(trim(nudge_flds(ifld))//' is not a 4d field in nudge_ncfile')
      if (dimids(1)/=dimid_lon .or. dimids(2)/=dimid_lat .or. dimids(4)/=dimid_time .or. &
        (dimids(3)/=dimid_lev .and. dimids(3)/=dimid_ilev)) &
        call shutdown(trim(nudge_flds(ifld))//' is not a valid 4d field in nudge_ncfile')

      if (dimids(3) == dimid_lev) then
        vcoord_is_midpoint(ifld) = .true.
      else
        vcoord_is_midpoint(ifld) = .false.
      endif
    enddo

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at closing '//trim(filename))

! if there are multiple external data files, the model will read the correct data file based on time
    nfile = count(len_trim(nudge_ncfile) > 0)
    allocate(nt(nfile))

    do ifile = 1, nfile
      filename = trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost)

      stat = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at opening '//trim(filename))

      stat = nf90_inq_dimid(ncid, 'time', dimid_time)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension time ID')

      stat = nf90_inquire_dimension(ncid, dimid_time, len=nt(ifile))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring dimension time length')

      if (nt(ifile) > maxnt) call shutdown('nudge_ncfile dimension time exceeds the maximum number of records')

      if (nudge_use_refdate) then
        stat = nf90_inq_varid(ncid, 'time', varid_time)
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable time ID')

        stat = nf90_get_var(ncid, varid_time, t(ifile, 1: nt(ifile)))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable time')
      else
        stat = nf90_inq_varid(ncid, 'date', varid_date)
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable date ID')

        stat = nf90_get_var(ncid, varid_date, dt(ifile, 1: nt(ifile)))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable date')

        stat = nf90_inq_varid(ncid, 'datesec', varid_datesec)
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at inquiring variable datesec ID')

        stat = nf90_get_var(ncid, varid_datesec, dts(ifile, 1: nt(ifile)))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable datesec')
      endif

      stat = nf90_close(ncid)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at closing '//trim(filename))
    enddo

    allocate(t0(nfile))
    allocate(t1(nfile))
    t0(1) = 1
    t1(1) = t0(1) + nt(1) - 1
    do ifile = 2, nfile
      t0(ifile) = t1(ifile-1) + 1
      t1(ifile) = t0(ifile) + nt(ifile) - 1
    enddo

    ntime = t1(nfile)
    if (ntime <= 1) call shutdown('nudge_ncfile must have at least 2 time points')

    allocate(time(ntime))
    start_datenum = to_datenum(start_year, start_day)
    if (nudge_use_refdate) then
      dt0 = (to_datenum(nudge_refdate(1), nudge_refdate(2)) - start_datenum) * 86400
      do ifile = 1, nfile
        time(t0(ifile): t1(ifile)) = t(ifile, 1: nt(ifile)) + dt0
      enddo
    else
      do ifile = 1, nfile
        do it = 1, nt(ifile)
          yr = dt(ifile, it) / 10000
          mn = mod(dt(ifile, it), 10000) / 100
          dy = mod(dt(ifile, it), 100)
          time(t0(ifile) + it - 1) = (to_datenum(yr, to_doy(yr, mn, dy)) - start_datenum) * 86400 + dts(ifile, it)
        enddo
      enddo
    endif

  end subroutine check
!-----------------------------------------------------------------------
  subroutine init
! synchronize external field paramters among processes and init subdomain grids

    use params_module, only: nlevp1, zibot, zpmid, zpint, glon0, glat, ispval
    use cons_module, only: dtr
    use input_module, only: mxlen_filename, nudge_ncpre, nudge_ncfile, nudge_ncpost, &
      nudge_flds, nudge_sponge, nudge_delta, nudge_power
    use fields_module, only: f4d_int=>f4d
    use mpi_module, only: mytid, lon0, lon1, lat0, lat1, TIEGCM_WORLD
    use netcdf, only: nf90_open, nf90_inq_varid, nf90_strerror, nf90_nowrite, nf90_noerr
    use mpi

    logical :: intersect
    integer :: ierror, stat, if4d, j, i, ls, rs
    real :: lb, rb, tb, bb, slon, latpart1, latpart2, lonpart
    character(len=mxlen_filename) :: filename
    integer, dimension(1) :: idx
    integer, dimension(13) :: request
    real, dimension(nlevp1) :: zout
    real, dimension(lon0: lon1) :: xc
    real, dimension(lat0: lat1) :: yc
    real, dimension(lon0: lon1, lat0: lat1) :: dist
    external :: shutdown

    call mpi_ibcast(nlon, 1, mpi_integer, 0, TIEGCM_WORLD, request(1), ierror)
    call mpi_ibcast(nlat, 1, mpi_integer, 0, TIEGCM_WORLD, request(2), ierror)
    call mpi_ibcast(nlev, 1, mpi_integer, 0, TIEGCM_WORLD, request(3), ierror)
    call mpi_ibcast(nilev, 1, mpi_integer, 0, TIEGCM_WORLD, request(4), ierror)
    call mpi_ibcast(ntime, 1, mpi_integer, 0, TIEGCM_WORLD, request(5), ierror)
    call mpi_ibcast(nfile, 1, mpi_integer, 0, TIEGCM_WORLD, request(6), ierror)
    call mpi_ibcast(nf4d, 1, mpi_integer, 0, TIEGCM_WORLD, request(7), ierror)
    call mpi_waitall(7, request(1: 7), mpi_statuses_ignore, ierror)

    if (mytid /= 0) then
      allocate(vcoord_is_midpoint(nf4d))
      allocate(time(ntime))
      allocate(nt(nfile))
      allocate(t0(nfile))
      allocate(t1(nfile))
      allocate(lon(nlon))
      allocate(lat(nlat))
      allocate(lev(nlev))
      allocate(ilev(nilev))
      allocate(f4d_idx(nf4d))
      allocate(varid(nf4d))
    endif

    call mpi_ibcast(wrap, 1, mpi_logical, 0, TIEGCM_WORLD, request(1), ierror)
    call mpi_ibcast(maxlev, 1, mpi_integer, 0, TIEGCM_WORLD, request(2), ierror)
    call mpi_ibcast(vcoord_is_midpoint, nf4d, mpi_logical, 0, TIEGCM_WORLD, request(3), ierror)
    call mpi_ibcast(time, ntime, mpi_integer, 0, TIEGCM_WORLD, request(4), ierror)
    call mpi_ibcast(nt, nfile, mpi_integer, 0, TIEGCM_WORLD, request(5), ierror)
    call mpi_ibcast(t0, nfile, mpi_integer, 0, TIEGCM_WORLD, request(6), ierror)
    call mpi_ibcast(t1, nfile, mpi_integer, 0, TIEGCM_WORLD, request(7), ierror)
    call mpi_ibcast(lon, nlon, mpi_real8, 0, TIEGCM_WORLD, request(8), ierror)
    call mpi_ibcast(lat, nlat, mpi_real8, 0, TIEGCM_WORLD, request(9), ierror)
    call mpi_ibcast(lev, nlev, mpi_real8, 0, TIEGCM_WORLD, request(10), ierror)
    call mpi_ibcast(ilev, nilev, mpi_real8, 0, TIEGCM_WORLD, request(11), ierror)
    call mpi_ibcast(lbc_flag, nlb, mpi_logical, 0, TIEGCM_WORLD, request(12), ierror)
    call mpi_ibcast(f4d_idx, nf4d, mpi_integer, 0, TIEGCM_WORLD, request(13), ierror)
    call mpi_waitall(13, request, mpi_statuses_ignore, ierror)

! smooth vertical transition from external fields to model fields, exponential decay
    allocate(vert_weight(maxlev, nf4d))
    do if4d = 1, nf4d
      if (trim(f4d_int(f4d_idx(if4d))%vcoord) == 'midpoints') then
        zout(1: maxlev) = zpmid(1: maxlev)
      else
        zout(1: maxlev) = zpint(1: maxlev)
      endif
      vert_weight(:, if4d) = exp(-((zout(1: maxlev) - zibot) / nudge_delta(2))**nudge_power(2))
    enddo

    do latbeg = lat0, lat1
      if (glat(latbeg) >= lat(1)) exit
    enddo
    do latend = lat1, lat0, -1
      if (glat(latend) <= lat(nlat)) exit
    enddo

    lonbeg = ispval
    lonend = ispval
    offbeg = ispval
    offend = ispval
    lonbeg1 = ispval
    lonend1 = ispval

    if (wrap) then
! if the external field covers the full longitude cycle, then the model subdomain is fully embedded
! move the edge of the external domain to cover the model subdomain
      slon = glon0((lon0+lon1)/2) - 180
      do offbeg = -1, 1
        if (lon(1)+offbeg*360<=slon .and. slon<=lon(nlon)+offbeg*360) exit
      enddo
      idx = minloc(abs(lon + offbeg*360 - slon))
      lonbeg = idx(1)
    else

! find the suitable 360 degree wrap of the external domain to intersect the model subdomain
! the following part finds the left edge of the model subdomain
      intersect = .false.
      left: do lonbeg = lon0, lon1
! the shift loop can be omitted if the external domain is also from -180 to 180 (shift=0)
! if the external domain is not from -180 to 180, the 360 degree wrap will be only one of -1,0,1
        do offbeg = -1, 1
          if (lon(1)+offbeg*360<=glon0(lonbeg) .and. glon0(lonbeg)<=lon(nlon)+offbeg*360) then
            intersect = .true.
            exit left
          endif
        enddo
      enddo left
      if (.not. intersect) then
        lonbeg = ispval
        offbeg = ispval
      endif

! similar for the right edge of the mode subdomain
      intersect = .false.
      right: do lonend = lon1, lon0, -1
        do offend = -1, 1
          if (lon(1)+offend*360<=glon0(lonend) .and. glon0(lonend)<=lon(nlon)+offend*360) then
            intersect = .true.
            exit right
          endif
        enddo
      enddo right
      if (.not. intersect) then
        lonend = ispval
        offend = ispval
      endif

! after the 360 degree wrap of the external domain, there will be three possibilities, discussed as follows:
! 1. lonbeg<lonend && offbeg==offend, left edge and right edge doesn't cross the 180 degree boundary,
! the overlapping region is (model_lon(lonbeg), model_lon(lonend)) and (lon(1)+offbeg*360, lon(nlon)+offend*360)
! 2. lonbeg<lonend && offbeg+1==offend, left edge doesn't cross the 180 degree boundary but right edge does,
! the overlapping region is splitted into two:
!   - 1. (model_lon(lonbeg), model_lon(lonend1)) vs (lon(1)+offbeg*360, lon(nlon)+offbeg*360)
!   - 2. (model_lon(lonbeg1), model_lon(lonend)) vs (lon(1)+offend*360, lon(nlon)+offend*360)
! the additional parameters lonbeg1,lonend1 of this situation are calculated below
! 3. lonbeg>lonend && offbeg==offend+1, left edge crosses the 180 degree boundary but right edge doesn't,
! left edge and right edge are swapped, the overlapping region also consists of two parts
!   - 1. (model_lon(lonbeg), model_lon(i1)) vs (lon(1)+offbeg*360, lon(nlon)+offbeg*360)
!   - 2. (model_lon(i0), model_lon(lonend)) vs (lon(1)+offend*360, lon(nlon)+offend*360)

      if (lonbeg<lonend .and. offbeg+1==offend) then
        do lonend1 = lonbeg, lonend
          if (glon0(lonend1) > lon(nlon)+offbeg*360) exit
        enddo
        lonend1 = lonend1 - 1
        do lonbeg1 = lonend, lonbeg, -1
          if (glon0(lonbeg1) < lon(1)+offend*360) exit
        enddo
        lonbeg1 = lonbeg1 + 1
      endif
    endif

! smooth horizontal transition from external fields to model fields
    if (latbeg < latend) then
      allocate(hori_weight(lon0: lon1, latbeg: latend))

! a sponge layer is imposed to allow smooth transition
      lb = lon(1) + nudge_sponge(1)
      rb = lon(nlon) - nudge_sponge(1)
      tb = lat(1) + nudge_sponge(1)
      bb = lat(nlat) - nudge_sponge(1)

! the sponge layer on latitudes
      do j = latbeg, latend
        if (glat(j) < tb) yc(j) = tb
        if (tb<=glat(j) .and. glat(j)<=bb) yc(j) = glat(j)
        if (glat(j) > bb) yc(j) = bb
      enddo

! the sponge layer on longitudes
      if (wrap) then
        xc = glon0(lon0: lon1)
      else
        do i = lon0, lon1
          do ls = -1, 1
            if (abs(lb+ls*360 - glon0(i)) <= 180) exit
          enddo
          do rs = -1, 1
            if (abs(rb+rs*360 - glon0(i)) <= 180) exit
          enddo
          if (glon0(i) < lb+ls*360) xc(i) = lb + ls*360
          if (lb+ls*360<=glon0(i) .and. glon0(i)<=rb+rs*360) xc(i) = glon0(i)
          if (glon0(i) > rb+rs*360) xc(i) = rb + rs*360
        enddo
      endif

! if the point is inside the inner square, the distance is zero
! if the point is outside the inner square, the distance is the great-circle distance to the nearest point
      do j = latbeg, latend
        latpart1 = sin((yc(j) - glat(j)) * dtr / 2)**2
        latpart2 = 1 - latpart1 - sin((yc(j) + glat(j)) * dtr / 2)**2
        do i = lon0, lon1
          lonpart = sin((xc(i) - glon0(i)) * dtr / 2)**2
          dist(i, j) = 2 * asin(sqrt(latpart1 + latpart2*lonpart))
        enddo
      enddo

! the relaxation function is an exponential function of the great-circle distance
      hori_weight = exp(-(dist(:, latbeg: latend) / (nudge_delta(1)*dtr))**nudge_power(1))
    endif

    ifile = 1
    itime = 1

    if (latbeg < latend) then
      allocate(lbc(lon0: lon1, latbeg: latend, 2, nlb))
      allocate(f4d(maxlev, lon0: lon1, latbeg: latend, 2, nf4d))

      filename = trim(nudge_ncpre)//trim(nudge_ncfile(1))//trim(nudge_ncpost)

      stat = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at opening '//trim(filename))

      do if4d = 1, nf4d
        stat = nf90_inq_varid(ncid, trim(nudge_flds(if4d)), varid(if4d))
        if (stat /= nf90_noerr) &
          call shutdown(trim(nf90_strerror(stat))//' at inquiring variable '//trim(nudge_flds(if4d))//' ID')
      enddo

      call read_data(itime, 1)
      call read_data(itime+1, 2)
    endif

  end subroutine init
!-----------------------------------------------------------------------
  subroutine update(modelsec)
! move the cursor itime to include modelsec inside the interval [time(itime), time(itime+1)]
! only processes with latbeg<latend calls this subroutine

    integer, intent(in) :: modelsec

    integer :: increment

    increment = 0
    do while (itime <= ntime-2)
      if (time(itime)<=modelsec .and. modelsec<=time(itime+1)) exit
      itime = itime + 1
      increment = increment + 1
    enddo

    if (increment >= 1) then
      if (increment == 1) then
! if itime is moved forward by 1, the previous time is then filled the the current time
        lbc(:, :, 1, :) = lbc(:, :, 2, :)
        f4d(:, :, :, 1, :) = f4d(:, :, :, 2, :)
      else
! only at the first step will both steps be re-calculated
        call read_data(itime, 1)
      endif
      call read_data(itime+1, 2)
    endif

  end subroutine update
!-----------------------------------------------------------------------
  subroutine read_data(itime, it)
! read the external field and interpolate to model grids
! only processes with latbeg<latend calls this subroutine

    use params_module, only: zibot, zpmid, zpint, glon0, glat
    use input_module, only: mxlen_filename, nudge_ncpre, nudge_ncfile, nudge_ncpost, nudge_flds, nudge_lbc, nudge_f4d
    use fields_module, only: f4d_int=>f4d
    use char_module, only: find_index
    use interp_module, only: interp3d
    use mpi_module, only: lon0, lon1
    use netcdf, only: nf90_close, nf90_open, nf90_inq_varid, nf90_get_var, nf90_strerror, nf90_nowrite, nf90_noerr

    integer, intent(in) :: itime, it

    integer :: stat, ifld, nk, ik, lb_idx
    character(len=mxlen_filename) :: filename
    real, dimension(nlon+1) :: lonp1
    real, dimension(max(nlev, nilev)) :: zin
    real, dimension(maxlev) :: zout
    real, dimension(1, lon0: lon1, latbeg: latend) :: lbc0
    real(kind=4), dimension(nlon, nlat, max(nlev, nilev)) :: f4d0
    real, dimension(max(nlev, nilev), nlon+1, nlat) :: ncf4d
    external :: shutdown

! open a new file to read if the current modelsec is over the last time index of the current file
    if (itime > t1(ifile)) then
      filename = trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost)

      stat = nf90_close(ncid)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at closing '//trim(filename))

      ifile = ifile + 1
      filename = trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost)

      stat = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at opening '//trim(filename))

      do ifld = 1, nf4d
        stat = nf90_inq_varid(ncid, trim(nudge_flds(ifld)), varid(ifld))
        if (stat /= nf90_noerr) &
          call shutdown(trim(nf90_strerror(stat))//' at inquiring variable '//trim(nudge_flds(ifld))//' ID')
      enddo
    endif

    if (wrap) then
      lonp1(1: nlon-lonbeg) = lon(lonbeg+1: nlon)
      lonp1(nlon-lonbeg+1: nlon+1) = lon(1: lonbeg+1) + 360
    endif

    do ifld = 1, nf4d
      if (vcoord_is_midpoint(ifld)) then
        nk = nlev
        zin(1: nk) = lev(1: nk)
      else
        nk = nilev
        zin(1: nk) = ilev(1: nk)
      endif

      stat = nf90_get_var(ncid, varid(ifld), f4d0(:, :, 1: nk), start=(/1, 1, 1, itime-t0(ifile)+1/), count=(/nlon, nlat, nk, 1/))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at getting variable '//trim(nudge_flds(ifld)))

      if (wrap) then
        do ik = 1, max(nlev, nilev)
          ncf4d(ik, 1: nlon-lonbeg, :) = f4d0(lonbeg+1: nlon, :, ik)
          ncf4d(ik, nlon-lonbeg+1: nlon, :) = f4d0(1: lonbeg, :, ik)
        enddo
        ncf4d(:, nlon+1, :) = ncf4d(:, 1, :)
      else
        do ik = 1, max(nlev, nilev)
          ncf4d(ik, 1: nlon, :) = f4d0(:, :, ik)
        enddo
      endif

! interpolate based on the subdomain intersection discussed in init
      if (nudge_lbc(ifld)) then
        lb_idx = find_index(nudge_flds(ifld), lb)
        if (lb_idx /= 0) then
          if (wrap) then
            lbc0 = interp3d((/zibot/), glon0(lon0: lon1), glat(latbeg: latend), &
              zin(1: nk), lonp1+offbeg*360, lat, ncf4d(1: nk, :, :))
            lbc(:, :, it, lb_idx) = lbc0(1, :, :)
          else
            if (lonbeg<lonend .and. offbeg==offend) then
              lbc0(:, lonbeg: lonend, :) = &
                interp3d((/zibot/), glon0(lonbeg: lonend), glat(latbeg: latend), &
                zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))
              lbc(lonbeg: lonend, :, it, lb_idx) = lbc0(1, lonbeg: lonend, :)
            endif

            if (lonbeg<lonend .and. offbeg+1==offend) then
              lbc0(:, lonbeg: lonend1, :) = &
                interp3d((/zibot/), glon0(lonbeg: lonend1), glat(latbeg: latend), &
                zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))
              lbc(lonbeg: lonend1, :, it, lb_idx) = lbc0(1, lonbeg: lonend1, :)
              lbc0(:, lonbeg1: lonend, :) = &
                interp3d((/zibot/), glon0(lonbeg1: lonend), glat(latbeg: latend), &
                zin(1: nk), lon+offend*360, lat, ncf4d(1: nk, 1: nlon, :))
              lbc(lonbeg1: lonend, :, it, lb_idx) = lbc0(1, lonbeg1: lonend, :)
            endif

            if (lonbeg>lonend .and. offbeg==offend+1) then
              lbc0(:, lonbeg: lon1, :) = &
                interp3d((/zibot/), glon0(lonbeg: lon1), glat(latbeg: latend), &
                zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))
              lbc(lonbeg: lon1, :, it, lb_idx) = lbc0(1, lonbeg: lon1, :)
              lbc0(:, lon0: lonend, :) = &
                interp3d((/zibot/), glon0(lon0: lonend), glat(latbeg: latend), &
                zin(1: nk), lon+offend*360, lat, ncf4d(1: nk, 1: nlon, :))
              lbc(lon0: lonend, :, it, lb_idx) = lbc0(1, lon0: lonend, :)
            endif
          endif
        endif
      endif

      if (nudge_f4d(ifld)) then
        if (trim(f4d_int(f4d_idx(ifld))%vcoord) == 'midpoints') then
          zout = zpmid(1: maxlev)
        else
          zout = zpint(1: maxlev)
        endif

        if (wrap) then
          f4d(:, :, :, it, ifld) = &
            interp3d(zout, glon0(lon0: lon1), glat(latbeg: latend), &
            zin(1: nk), lonp1+offbeg*360, lat, ncf4d(1: nk, :, :))
        else
          if (lonbeg<lonend .and. offbeg==offend) &
            f4d(:, lonbeg: lonend, :, it, ifld) = &
            interp3d(zout, glon0(lonbeg: lonend), glat(latbeg: latend), &
            zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))

          if (lonbeg<lonend .and. offbeg+1==offend) then
            f4d(:, lonbeg: lonend1, :, it, ifld) = &
              interp3d(zout, glon0(lonbeg: lonend1), glat(latbeg: latend), &
              zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))
            f4d(:, lonbeg1: lonend, :, it, ifld) = &
              interp3d(zout, glon0(lonbeg1: lonend), glat(latbeg: latend), &
              zin(1: nk), lon+offend*360, lat, ncf4d(1: nk, 1: nlon, :))
          endif

          if (lonbeg>lonend .and. offbeg==offend+1) then
            f4d(:, lonbeg: lon1, :, it, ifld) = &
              interp3d(zout, glon0(lonbeg: lon1), glat(latbeg: latend), &
              zin(1: nk), lon+offbeg*360, lat, ncf4d(1: nk, 1: nlon, :))
            f4d(:, lon0: lonend, :, it, ifld) = &
              interp3d(zout, glon0(lon0: lonend), glat(latbeg: latend), &
              zin(1: nk), lon+offend*360, lat, ncf4d(1: nk, 1: nlon, :))
          endif
        endif
      endif
    enddo

  end subroutine read_data
!-----------------------------------------------------------------------
  subroutine finalize
! only processes with latbeg<latend calls this subroutine

    use netcdf, only: nf90_close, nf90_strerror, nf90_noerr

    integer :: stat
    external :: shutdown

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at closing nudge_ncfile')

  end subroutine finalize
!-----------------------------------------------------------------------
end module nudge_module
