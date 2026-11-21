import os

from numpy import (ndarray, interp, log, exp, linspace, allclose, mean, isfinite, maximum, errstate,
                   asarray, count_nonzero)
from netCDF4 import Dataset

def interp2d(variable, inlat, inlon, outlat, outlon):
    """Interpolate a (lat, lon) field to the output grid; longitude is periodic."""

    ninlat = len(inlat)
    noutlon = len(outlon)

    var0 = ndarray(shape=(ninlat, noutlon))
    for ilat in range(ninlat):
        var0[ilat, :] = interp(x=outlon, xp=inlon, fp=variable[ilat, :], period=360)

    var1 = ndarray(shape=(len(outlat), noutlon))
    for ilon in range(noutlon):
        var1[:, ilon] = interp(x=outlat, xp=inlat, fp=var0[:, ilon])

    return var1


def interp3d(variable, inlev, inlat, inlon, outlev, outlat, outlon, extrap):
    """Interpolate a (lev, lat, lon) field; `extrap` above the top input level is 'constant',
    'linear' or 'exponential'."""
    
    ninlev = len(inlev)
    noutlev = len(outlev)
    noutlat = len(outlat)
    noutlon = len(outlon)

    var0 = ndarray(shape=(ninlev, noutlat, noutlon))
    for ik in range(ninlev):
        var0[ik, :, :] = interp2d(variable=variable[ik, :, :], inlat=inlat, inlon=inlon, outlat=outlat, outlon=outlon)

    # lastidx: number of output levels within the input range
    for lastidx in range(noutlev):
        if outlev[lastidx] > inlev[-1]:
            break

    if lastidx == noutlev-1 and outlev[lastidx] <= inlev[-1]:
        lastidx = noutlev

    v1 = ndarray(shape=noutlev)
    var1 = ndarray(shape=(noutlev, noutlat, noutlon))
    for ilat in range(noutlat):
        for ilon in range(noutlon):
            v0 = var0[:, ilat, ilon]
            if extrap == 'constant':
                v1[0: lastidx] = interp(x=outlev[0: lastidx], xp=inlev, fp=v0)
                if lastidx < noutlev:
                    v1[lastidx: noutlev] = v1[lastidx-1]
            elif extrap == 'linear':
                v1 = interp_linear(v0, inlev, outlev, lastidx)
            elif extrap == 'exponential':
                # a column with negative values (e.g. N2D near the top) has no log: fall back to
                # linear, floored at 0
                with errstate(divide='ignore', invalid='ignore', over='ignore'):
                    v1 = exp(interp_linear(log(v0), inlev, outlev, lastidx))
                if not isfinite(v1).all():
                    v1 = maximum(interp_linear(v0, inlev, outlev, lastidx), 0)
            else:
                exit('Extrapolation method must be one of constant/linear/exponential')
            var1[:, ilat, ilon] = v1

    return var1


def interp_linear(v0, inlev, outlev, lastidx):
    """Interpolate a column linearly; levels from `lastidx` up are extrapolated from the top two."""
    ninlev = len(inlev)
    noutlev = len(outlev)
    v1 = ndarray(shape=noutlev)
    v1[0: lastidx] = interp(x=outlev[0: lastidx], xp=inlev, fp=v0)
    k = (v0[ninlev-1] - v0[ninlev-2]) / (inlev[ninlev-1] - inlev[ninlev-2])
    if lastidx < noutlev:
        v1[lastidx: noutlev] = k * (outlev[lastidx: noutlev] - inlev[ninlev-1]) + v0[ninlev-1]
    return v1


def interpic(fin, hres, vres, zitop, fout):
    """Regrid TIE-GCM primary history `fin` to a new grid and write it to `fout`."""
    print(f"Interpolating primary file {fin} to create new primary file {fout} at horizontal resolution {hres} and vertical resolution {vres} with zitop {zitop}.")
    lower_cap = 1e-12
    fill_top = ['TN', 'UN', 'VN', 'OP', 'TI', 'TE', 'N2D', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'OP_NM', 'PED', 'HALL']
    mixing_ratio = ['O2', 'O1', 'HE', 'N2D', 'N4S', 'NO', 'AR', 'O2_NM', 'O1_NM', 'HE_NM', 'N2D_NM', 'N4S_NM', 'NO_NM', 'AR_NM']
    extrap_method = {'TN': 'exponential', 'UN': 'linear', 'VN': 'linear',
            'O2': 'exponential', 'O1': 'exponential', 'HE': 'exponential',
            'OP': 'exponential', 'N2D': 'exponential', 'N4S': 'exponential',
            'NO': 'exponential', 'AR': 'exponential', 'TI': 'exponential',
            'TE': 'exponential', 'NE': 'exponential', 'OMEGA': 'linear', 
            'O2P': 'constant', 'Z': 'exponential', 'POTEN': 'linear',
            'TN_NM': 'exponential', 'UN_NM': 'linear', 'VN_NM': 'linear',
            'O2_NM': 'exponential', 'O1_NM': 'exponential', 'HE_NM': 'exponential',
            'OP_NM': 'exponential', 'N2D_NM': 'exponential', 'N4S_NM': 'exponential',
            'NO_NM': 'exponential', 'AR_NM': 'exponential', 'PED': 'exponential', 'HALL': 'exponential'}

    nlon = int(360 / hres)
    lon = linspace(start=-180, stop=180-hres, num=nlon)
    nlat = int(180 / hres)
    lat = linspace(start=-90+hres/2, stop=90-hres/2, num=nlat)
    nlev = int((zitop + 7) / vres) + 1
    ilev = linspace(start=-7, stop=zitop, num=nlev)
    lev = ilev + vres/2

    src = Dataset(filename=fin)
    dst = Dataset(filename=fout, mode='w')

    print('Creating new primary file: ', fout)

    for name in src.ncattrs():
        setattr(dst, name, getattr(src, name))

    for dimname, dimension in src.dimensions.items():
        if dimname == 'time':
            nt = dimension.size
            dst.createDimension(dimname='time')
        elif dimname == 'lon':
            dst.createDimension(dimname='lon', size=nlon)
        elif dimname == 'lat':
            dst.createDimension(dimname='lat', size=nlat)
        elif dimname == 'lev':
            dst.createDimension(dimname='lev', size=nlev)
        elif dimname == 'ilev':
            dst.createDimension(dimname='ilev', size=nlev)
        elif dimname == 'mtimedim':
            dst.createDimension(dimname='mtimedim', size=4)
        else:
            dst.createDimension(dimname=dimname, size=dimension.size)

    lon0 = src['lon'][:]
    lat0 = src['lat'][:]
    lev0 = src['lev'][:]
    ilev0 = src['ilev'][:]

    nlon0 = len(lon0)
    nlat0 = len(lat0)
    nlev0 = len(lev0)

    # add the poles: the boundary latitudes shift with horizontal resolution
    lat0_bnd = ndarray(shape=nlat0+2)
    lat0_bnd[0] = -90
    lat0_bnd[1: nlat0+1] = lat0
    lat0_bnd[nlat0+1] = 90

    for varname, variable in src.variables.items():
        if varname == 'coupled_cmit':
            varout = dst.createVariable(varname='coupled_mage', datatype=variable.datatype, dimensions=variable.dimensions)
        else:
            varout = dst.createVariable(varname=varname, datatype=variable.datatype, dimensions=variable.dimensions)

        for name in variable.ncattrs():
            setattr(varout, name, getattr(variable, name))

        if varname == 'time':
            varout[:] = variable[:]
        elif varname == 'lon':
            varout[:] = lon
        elif varname == 'lat':
            varout[:] = lat
        elif varname == 'lev':
            varout[:] = lev
        elif varname == 'ilev':
            varout[:] = ilev

        # only in some non-standard histories
        elif varname in ['lat_bnds', 'lon_bnds', 'gw', 'area']:
            continue

        # old 3-element mtime -> 4
        elif varname == 'mtime':
            if src.dimensions['mtimedim'].size == 3:
                varout[:, 0: 3] = variable[:]
                varout[:, 3] = 0
            else:
                varout[:] = variable[:]

        # source run had calc_helium == 0
        elif varname in ['HE', 'HE_NM'] and allclose(variable, 0):
            varout[:] = lower_cap

        elif variable.dimensions == ('time', 'lat', 'lon'):
            var2d_bnd = ndarray(shape=(nlat0+2, nlon0))
            for it in range(nt):
                # pole value = mean of the highest latitude circle
                var2d_bnd[0, :] = mean(variable[it, 0, :])
                var2d_bnd[1: nlat0+1, :] = variable[it, :, :]
                var2d_bnd[nlat0+1, :] = mean(variable[it, nlat0-1, :])
                varout[it, :, :] = interp2d(variable=var2d_bnd, inlat=lat0_bnd, inlon=lon0, outlat=lat, outlon=lon)

        elif len(variable.dimensions) == 4:
            if not variable.dimensions in [('time', 'lev', 'lat', 'lon'), ('time', 'ilev', 'lat', 'lon')]:
                exit('Invalid 4d field: '+varname)

            if variable.dimensions[1] == 'lev':
                levin = lev0
                levout = lev
            else:
                levin = ilev0
                levout = ilev

            # top level holds fill values for these
            if varname in fill_top:
                nlevin = nlev0 - 1
            else:
                nlevin = nlev0

            var3d_bnd = ndarray(shape=(nlevin, nlat0+2, nlon0))
            for it in range(nt):
                for ik in range(nlevin):
                    var3d_bnd[ik, 0, :] = mean(variable[it, ik, 0, :])
                    var3d_bnd[ik, 1: nlat0+1, :] = variable[it, ik, :, :]
                    var3d_bnd[ik, nlat0+1, :] = mean(variable[it, ik, nlat0-1, :])
                varout[it, :, :, :] = interp3d(variable=var3d_bnd, inlev=levin[0: nlevin], inlat=lat0_bnd, inlon=lon0,
                    outlev=levout, outlat=lat, outlon=lon, extrap=extrap_method.get(varname, 'constant'))

            # mixing ratios in [0, 1], major species floored at lower_cap
            if varname in mixing_ratio:
                v = varout[:]
                lower = lower_cap if varname in ['O2', 'O1', 'HE', 'O2_NM', 'O1_NM', 'HE_NM'] else 0
                v[v < lower] = lower
                v[v > 1] = 1
                varout[:] = v

        else:
            varout[:] = variable[:]

    # renormalize the major species with an interpolated N2
    for ext in ['', '_NM']:
        if 'HE'+ext in src.variables.keys():
            N2 = 1 - src['O2'+ext][:] - src['O1'+ext][:] - src['HE'+ext][:]
        else:
            N2 = 1 - src['O2'+ext][:] - src['O1'+ext][:]
            dst.createVariable(varname='HE'+ext, datatype='f8', dimensions=('time', 'lev', 'lat', 'lon'))
            dst['HE'+ext][:] = lower_cap

        N2n = ndarray(shape=(nt, nlev, nlat, nlon))
        N2_bnd = ndarray(shape=(nlev0, nlat0+2, nlon0))
        for it in range(nt):
            for ik in range(nlev0):
                N2_bnd[ik, 0, :] = mean(N2[it, ik, 0, :])
                N2_bnd[ik, 1: nlat0+1, :] = N2[it, ik, :, :]
                N2_bnd[ik, nlat0+1, :] = mean(N2[it, ik, nlat0-1, :])
            N2n[it, :, :, :] = interp3d(variable=N2_bnd, inlev=lev0, inlat=lat0_bnd, inlon=lon0,
                outlev=lev, outlat=lat, outlon=lon, extrap='exponential')

        N2n[N2n < lower_cap] = lower_cap
        N2n[N2n > 1] = 1
        O2n = dst['O2'+ext][:]
        O1n = dst['O1'+ext][:]
        HEn = dst['HE'+ext][:]

        normalize = O2n + O1n + HEn + N2n
        dst['O2'+ext][:] = O2n / normalize
        dst['O1'+ext][:] = O1n / normalize
        dst['HE'+ext][:] = HEn / normalize

    # New 2D variables since TIEGCM v3.0
    for varname in ['gzigm1', 'gzigm2', 'gnsrhs']:
        if not varname in src.variables.keys():
            newvarout = dst.createVariable(varname=varname, datatype='f8', dimensions=('time', 'lat', 'lon'))
            newvarout[:] = 0

    # the model only reports NaNs in a source history and runs on, so refuse them here
    nonfinite = {}
    for varname, variable in dst.variables.items():
        if variable.dtype.kind == 'f':
            n = count_nonzero(~isfinite(asarray(variable[:])))
            if n:
                nonfinite[varname] = int(n)

    src.close()
    dst.close()
    if nonfinite:
        os.remove(fout)
        raise ValueError(f"interpic: the regrid of {fin} has non-finite values (field: count) {nonfinite}; "
                         f"{fout} was not written. Check the source file for NaN/Inf.")
