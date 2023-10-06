#!/usr/bin/env python
from numpy import ndarray, interp, log, exp, linspace, allclose, maximum
from sys import exit
from netCDF4 import Dataset
import argparse
from argparse import RawTextHelpFormatter

# interpolate old or new TIEGCM primary history file into v3.0 format
# support any combinations of horizontal/vertical resolutions and upper boundaries


def interp2d(variable, inlat, inlon, outlat, outlon):
    ninlat = len(inlat)
    noutlon = len(outlon)

    var0 = ndarray(shape=(ninlat, noutlon))
    for i in range(ninlat):
        var0[i, :] = interp(x=outlon, xp=inlon, fp=variable[i, :], period=360)

    var1 = ndarray(shape=(len(outlat), noutlon))
    for i in range(noutlon):
        var1[:, i] = interp(x=outlat, xp=inlat, fp=var0[:, i])

    return var1


def interp3d(variable, inlev, inlat, inlon, outlev, outlat, outlon, extrap):
    ninlev = len(inlev)
    noutlev = len(outlev)
    noutlat = len(outlat)
    noutlon = len(outlon)

    var0 = ndarray(shape=(ninlev, noutlat, noutlon))
    for i in range(ninlev):
        var0[i, :, :] = interp2d(variable=variable[i, :, :], inlat=inlat, inlon=inlon, outlat=outlat, outlon=outlon)

    for lastidx in range(noutlev):
        if outlev[lastidx] > inlev[-1]:
            break

    var1 = ndarray(shape=(noutlev, noutlat, noutlon))
    if extrap == 'constant':
        for i in range(noutlat):
            for j in range(noutlon):
                var1[:, i, j] = interp(x=outlev, xp=inlev, fp=var0[:, i, j])
    elif extrap == 'linear':
        for i in range(noutlat):
            for j in range(noutlon):
                var1[0: lastidx, i, j] = interp(x=outlev[0: lastidx], xp=inlev, fp=var0[:, i, j])
                k = (var0[ninlev-1, i, j] - var0[ninlev-2, i, j]) / (inlev[ninlev-1] - inlev[ninlev-2])
                var1[lastidx: noutlev, i, j] = k * (outlev[lastidx: noutlev] - inlev[ninlev-1]) + var0[ninlev-1, i, j]
    elif extrap == 'exponential':
        v1 = ndarray(shape=noutlev)
        for i in range(noutlat):
            for j in range(noutlon):
                v0 = log(var0[:, i, j])
                v1[0: lastidx] = interp(x=outlev[0: lastidx], xp=inlev, fp=v0)
                k = (v0[ninlev-1] - v0[ninlev-2]) / (inlev[ninlev-1] - inlev[ninlev-2])
                v1[lastidx: noutlev] = k * (outlev[lastidx: noutlev] - inlev[ninlev-1]) + v0[ninlev-1]
                var1[:, i, j] = exp(v1)
    else:
        exit('Extrapolation method must be one of constant/linear/exponential')

    return var1

if __name__ == "__main__":

    fnIn = None
    hori = 2.5
    vert = 0.25
    zitop = 5
    MainS = """ This script will scale a TIEGCM primary restart file to a new compatible TIEGCM v3.0 restart file.
        """
    
    #Parse commandline input
    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-hres',type=float,metavar="HoriRes",default=hori,help="Horizontal Resolution (default: %(default)s)")
    parser.add_argument('-vres',type=float,metavar="VertRes",default=vert,help="Vertical Resolution (default: %(default)s)")
    parser.add_argument('-zitop',type=int,metavar="ZPMax",default=zitop,help="Top Pressure Level (default: %(default)s)")
    parser.add_argument('-fin', type=str,default=None,help="Required: Input NetCDF primary file to be scaled")
    parser.add_argument('-fout', type=str,default=None,help="Filename for Scaled NetCDF primary file (default: {fIn}_{hres}x{vres}_{zitop}.nc)")
    
    #Finalize parsing
    args = parser.parse_args()
    
    #Check for required file
    if (args.fin is None): raise Exception('Error: -fin file needs to be specified')
    
    fnIn  = args.fin
    hori  = args.hres
    vert  = args.vres
    zitop = args.zitop
    fnOut = args.fout

    #Do Stuff
    fill_top = ['TN', 'UN', 'VN', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'OP_NM']
    mixing_ratio = ['O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'N2D', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'HE_NM', 'AR_NM']
    extrap_method = {'TN': 'exponential', 'UN': 'linear', 'VN': 'linear', 'O2': 'exponential', 'O1': 'exponential', 'N4S': 'exponential',
            'NO': 'exponential', 'HE': 'exponential', 'AR': 'exponential', 'OP': 'exponential', 'N2D': 'constant', 'TI': 'exponential',
            'TE': 'exponential', 'NE': 'exponential', 'O2P': 'constant', 'OMEGA': 'linear', 'Z': 'exponential', 'POTEN': 'linear',
            'TN_NM': 'exponential', 'UN_NM': 'linear', 'VN_NM': 'linear', 'O2_NM': 'exponential', 'O1_NM': 'exponential',
            'N4S_NM': 'exponential', 'NO_NM': 'exponential', 'OP_NM': 'exponential', 'HE_NM': 'exponential', 'AR_NM': 'exponential'}

    nlon = int(360 / hori)
    lon = linspace(start=-180, stop=180-hori, num=nlon)
    nlat = int(180 / hori)
    lat = linspace(start=-90+hori/2, stop=90-hori/2, num=nlat)
    nlev = int((zitop + 7) / vert) + 1
    ilev = linspace(start=-7, stop=zitop, num=nlev)
    lev = ilev + vert/2

    src = Dataset(filename=fnIn)
    
    if (fnOut is None):
        dst = Dataset(filename=fnIn+'_{:g}x{:g}_z{:g}.nc'.format(hori, vert, zitop), mode='w')
    else:
        if (fnOut[-2:-1] != 'nc'): fnOut = fnOut + '.nc'
        dst = Dataset(filename=fnOut,mode='w')

    print("Creating new primary file: ",fnOut)

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

    print("Copying and Scaling variables now. Grab some coffee while you wait.")

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

        elif varname == 'mtime':
            if src.dimensions['mtimedim'].size == 3:
                varout[:, 0: 3] = variable[:]
                varout[:, 3] = 0
            else:
                varout[:] = variable[:]

        elif varname in ['HE', 'HE_NM'] and allclose(variable, 0):
            varout[:] = 0

        elif variable.dimensions == ('time', 'lat', 'lon'):
            for i in range(nt):
                varout[i, :, :] = interp2d(variable=variable[i, :, :], inlat=lat0, inlon=lon0, outlat=lat, outlon=lon)

        elif len(variable.dimensions) == 4:
            if variable.dimensions == ('time', 'lev', 'lat', 'lon'):
                if varname in fill_top:
                    levin = lev0[:-1]
                    v = variable[:, :-1, :, :]
                else:
                    levin = lev0
                    v = variable
                levout = lev
            elif variable.dimensions == ('time', 'ilev', 'lat', 'lon'):
                if varname in fill_top:
                    levin = ilev0[:-1]
                    v = variable[:, :-1, :, :]
                else:
                    levin = ilev0
                    v = variable
                levout = ilev
            else:
                exit('Invalid 4d field: '+varname)

            for i in range(nt):
                varout[i, :, :, :] = interp3d(variable=v[i, :, :, :], inlev=levin, inlat=lat0, inlon=lon0, outlev=levout, outlat=lat, outlon=lon, extrap=extrap_method[varname])

            if varname in mixing_ratio:
                v = varout[:]
                v[v > 1] = 1
                varout[:] = v

        else:
            varout[:] = variable[:]

    for ext in ['', '_NM']:
        try:
            N2 = maximum(1 - src['O2'+ext][:] - src['O1'+ext][:] - src['HE'+ext][:], 1e-6)
        except:
            N2 = maximum(1 - src['O2'+ext][:] - src['O1'+ext][:] , 1e-6)

        N2n = ndarray(shape=(nt, nlev, nlat, nlon))
        for i in range(nt):
            N2n[i, :, :, :] = interp3d(variable=N2[i, :, :, :], inlev=lev0, inlat=lat0, inlon=lon0, outlev=lev, outlat=lat, outlon=lon, extrap='exponential')
        O1n = dst['O1'+ext][:]
        try:
            HEn = dst['HE'+ext][:]
        except:
            HEn = 0.
        scale = (1 - dst['O2'+ext][:] - N2n) / (O1n + HEn)
        dst['O1'+ext][:] = O1n * scale
        try:
            dst['HE'+ext][:] = HEn * scale
        except:
            dst.createVariable(varname='HE'+ext, datatype='f8', dimensions=('time','lev','lat','lon'))
            dst['HE'+ext][:] = 0.

    for varname in ['gzigm1', 'gzigm2', 'gnsrhs']:
        if not varname in src.variables.keys():
            newvarout = dst.createVariable(varname=varname, datatype='f8', dimensions=('time', 'lat', 'lon'))
            newvarout[:] = 0

    src.close()
    dst.close()
    
    print("We're done now, so closing out the files")
