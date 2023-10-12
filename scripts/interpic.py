#!/usr/bin/env python

from numpy import ndarray, interp, log, exp, linspace, allclose, mean
from sys import exit
from netCDF4 import Dataset
from argparse import ArgumentParser
from os.path import isfile, splitext


# interpolate old or new TIEGCM primary history file into v3.0 format
# support any combinations of horizontal/vertical resolutions and upper boundary levels


def interp2d(variable, inlat, inlon, outlat, outlon):
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
    ninlev = len(inlev)
    noutlev = len(outlev)
    noutlat = len(outlat)
    noutlon = len(outlon)

    var0 = ndarray(shape=(ninlev, noutlat, noutlon))
    for ik in range(ninlev):
        var0[ik, :, :] = interp2d(variable=variable[ik, :, :], inlat=inlat, inlon=inlon, outlat=outlat, outlon=outlon)

    # Find the last index of outlev falling in the range of inlev
    for lastidx in range(noutlev):
        if outlev[lastidx] > inlev[-1]:
            break

    # If outlev is completely embedded in inlev (interpolation only), the end point needs to be added separately
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
                v1[0: lastidx] = interp(x=outlev[0: lastidx], xp=inlev, fp=v0)
                k = (v0[ninlev-1] - v0[ninlev-2]) / (inlev[ninlev-1] - inlev[ninlev-2])
                if lastidx < noutlev:
                    v1[lastidx: noutlev] = k * (outlev[lastidx: noutlev] - inlev[ninlev-1]) + v0[ninlev-1]
            elif extrap == 'exponential':
                v0 = log(v0)
                v1[0: lastidx] = interp(x=outlev[0: lastidx], xp=inlev, fp=v0)
                k = (v0[ninlev-1] - v0[ninlev-2]) / (inlev[ninlev-1] - inlev[ninlev-2])
                if lastidx < noutlev:
                    v1[lastidx: noutlev] = k * (outlev[lastidx: noutlev] - inlev[ninlev-1]) + v0[ninlev-1]
                v1 = exp(v1)
            else:
                exit('Extrapolation method must be one of constant/linear/exponential')
            var1[:, ilat, ilon] = v1

    return var1


def interpic(fin, hres, vres, zitop, fout):
    # Some additional attributes for 4D fields
    lower_cap = 1e-6
    fill_top = ['TN', 'UN', 'VN', 'OP', 'TI', 'TE', 'N2D', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'OP_NM']
    mixing_ratio = ['O2', 'O1', 'HE', 'N4S', 'NO', 'AR', 'N2D', 'O2_NM', 'O1_NM', 'HE_NM', 'N4S_NM', 'NO_NM', 'AR_NM']
    extrap_method = {'TN': 'exponential', 'UN': 'linear', 'VN': 'linear', 'O2': 'exponential', 'O1': 'exponential', 'HE': 'exponential',
            'OP': 'exponential', 'N4S': 'exponential', 'NO': 'exponential', 'AR': 'exponential', 'TI': 'exponential', 'TE': 'exponential',
            'NE': 'exponential', 'OMEGA': 'linear', 'N2D': 'constant',  'O2P': 'constant', 'Z': 'exponential', 'POTEN': 'linear',
            'TN_NM': 'exponential', 'UN_NM': 'linear', 'VN_NM': 'linear', 'O2_NM': 'exponential', 'O1_NM': 'exponential', 'HE_NM': 'exponential',
            'OP_NM': 'exponential', 'N4S_NM': 'exponential', 'NO_NM': 'exponential', 'AR_NM': 'exponential'}

    nlon = int(360 / hres)
    lon = linspace(start=-180, stop=180-hres, num=nlon)
    nlat = int(180 / hres)
    lat = linspace(start=-90+hres/2, stop=90-hres/2, num=nlat)
    nlev = int((zitop + 7) / vres) + 1
    ilev = linspace(start=-7, stop=zitop, num=nlev)
    lev = ilev + vres/2

    src = Dataset(filename=fin)
    dst = Dataset(filename=fout, mode='w')

    print("Creating new primary file: ",fout)

    # Copy all attributes from old to new files (even though many of them are not actually used)
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

    # Bound latitudes with two poles since the change of horizontal resolutions lead to a boundary latitude shift
    lat0_bnd = ndarray(shape=nlat0+2)
    lat0_bnd[0] = -90
    lat0_bnd[1: nlat0+1] = lat0
    lat0_bnd[nlat0+1] = 90

    for varname, variable in src.variables.items():
        # Name change only
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

        # Change from old format (3 digits) to new format (4 digits)
        elif varname == 'mtime':
            if src.dimensions['mtimedim'].size == 3:
                varout[:, 0: 3] = variable[:]
                varout[:, 3] = 0
            else:
                varout[:] = variable[:]

        # If the old file was from a run with calc_helium==0, then set a constant for Helium in the new file (don't interpolate)
        elif varname in ['HE', 'HE_NM'] and allclose(variable, 0):
            varout[:] = lower_cap

        # 3D fields
        elif variable.dimensions == ('time', 'lat', 'lon'):
            var2d_bnd = ndarray(shape=(nlat0+2, nlon0))
            for it in range(nt):
                # Set pole fields as the average of the highest latitude circle
                var2d_bnd[0, :] = mean(variable[it, 0, :])
                var2d_bnd[1: nlat0+1, :] = variable[it, :, :]
                var2d_bnd[nlat0+1, :] = mean(variable[it, nlat0-1, :])
                varout[it, :, :] = interp2d(variable=var2d_bnd, inlat=lat0_bnd, inlon=lon0, outlat=lat, outlon=lon)

        # 4D fields
        elif len(variable.dimensions) == 4:
            if not variable.dimensions in [('time', 'lev', 'lat', 'lon'), ('time', 'ilev', 'lat', 'lon')]:
                exit('Invalid 4d field: '+varname)

            if variable.dimensions[1] == 'lev':
                levin = lev0
                levout = lev
            else:
                levin = ilev0
                levout = ilev

            # If the topmost level are filling values, exclude that level
            if varname in fill_top:
                nlevin = nlev0 - 1
            else:
                nlevin = nlev0

            var3d_bnd = ndarray(shape=(nlevin, nlat0+2, nlon0))
            for it in range(nt):
                # Set pole fields as the average of the highest latitude circle
                for ik in range(nlevin):
                    var3d_bnd[ik, 0, :] = mean(variable[it, ik, 0, :])
                    var3d_bnd[ik, 1: nlat0+1, :] = variable[it, ik, :, :]
                    var3d_bnd[ik, nlat0+1, :] = mean(variable[it, ik, nlat0-1, :])
                varout[it, :, :, :] = interp3d(variable=var3d_bnd, inlev=levin[0: nlevin], inlat=lat0_bnd, inlon=lon0,
                    outlev=levout, outlat=lat, outlon=lon, extrap=extrap_method[varname])

            # Mixing ratio must lie within [0, 1], note that the exponential extrapolation gurantees positivity
            if varname in mixing_ratio:
                v = varout[:]
                # In addition, major species have a lower cap
                if varname in ['O2', 'O1', 'HE', 'O2_NM', 'O1_NM', 'HE_NM']:
                    v[v < lower_cap] = lower_cap
                v[v > 1] = 1
                varout[:] = v

        else:
            varout[:] = variable[:]

    # N2 needs to be extrapolated to check the validity of other major species (O2, O1, HE)
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

    src.close()
    dst.close()


if __name__ == '__main__':
    # Commandline options
    parser = ArgumentParser(description='This script will scale a TIEGCM primary restart file to a new compatible TIEGCM v3.0 restart file.')
    parser.add_argument('-fin', type=str, required=True, help='Input NetCDF primary history to be scaled (Required)')
    parser.add_argument('-hres', type=float, help='Horizontal resolution (Default: Read from fin)')
    parser.add_argument('-vres', type=float, help='Vertical resolution (Default: Read from fin)')
    parser.add_argument('-zitop', type=float, help='Top pressure level (Default: Read from fin)')
    parser.add_argument('-fout', type=str, help='Filename for scaled NetCDF primary history (Default: Resolution appended to fin)')

    # Parse commandline input
    args = parser.parse_args()

    # Check for required file
    if (not isfile(args.fin)):
        exit('Error: Input file does not exist')

    # Obtain arguments from input file
    src = Dataset(filename=args.fin)
    if (args.hres is None):
        lat0 = src['lat'][0: 2]
        args.hres = lat0[1] - lat0[0]
    if (args.vres is None):
        lev0 = src['lev'][0: 2]
        args.vres = lev0[1] - lev0[0]
    if (args.zitop is None):
        args.zitop = src['ilev'][-1]
    src.close()

    # Setup default name for output file
    if (args.fout is None):
        args.fout = '{:s}_{:g}x{:g}_z{:g}.nc'.format(splitext(args.fin)[0], args.hres, args.vres, args.zitop)

    # Echo all information
    print('Input file: {:s}\nHorizontal resolution: {:g}\nVertical resolution: {:g}\nTop pressure level: {:g}\nOutput file: {:s}'
        .format(args.fin, args.hres, args.vres, args.zitop, args.fout))

    # Do interpolation
    print('Processing...')
    interpic(args.fin, args.hres, args.vres, args.zitop, args.fout)
    print('Finished')
