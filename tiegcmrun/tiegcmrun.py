#!/usr/bin/env python


"""makeitso for the TIE-GCM software.

This script is performs all of the steps that are required to prepare to run a TIE-GCM simulation run. By default, this script is interactive - the user
is prompted for each decision  that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BENCH" - The user requests for a benchmark run for TIE-GCM 

"BASIC" (the default) - the user is prompted to set only a small subset of MAGE
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.
"""

# Import standard modules
import argparse
from argparse import ArgumentParser
import copy
import datetime
from datetime import datetime, timedelta
import json
import os
from os.path import isfile, splitext
import subprocess
import sys
import shutil
import filecmp
import fnmatch
from textwrap import dedent
from fractions import Fraction
from math import ceil

# Import 3rd-party modules
#import numpy as np
from numpy import pad, ndarray, interp, log, exp, linspace, allclose, mean
from netCDF4 import Dataset
import xarray as xr
from jinja2 import Template



# Program constants
RED = '\033[31m'  # Red text
GREEN = '\033[32m'  # Green text
YELLOW = '\033[33m'  # Yellow text
RESET = '\033[0m'  # Reset to default color


# Program description.
DESCRIPTION = "Interactive script to prepare a TIEGCM model run."

# Indent level for JSON output.
JSON_INDENT = 4

# Path to current tiegcm datafiles
try:
    TIEGCMDATA = os.environ["TIEGCMDATA"]
except:
    os.environ['TIEGCMDATA'] = input(f'{RED}Unable to get TIEGCMDATA environment variable.{RESET}\n{YELLOW}Use command "export TIEGCMDATA=Path/To/Data" to set environment variable.{RESET}\nEnter TIEGCM data directory: ')
    TIEGCMDATA = os.environ["TIEGCMDATA"]

# Path to current tiegcm installation
try:
    TIEGCMHOME = os.environ["TIEGCMHOME"]
except:
    os.environ['TIEGCMHOME'] = input(f'{RED}Unable to get TIEGCMHOME environment variable.{RESET}\n{YELLOW}Use command "export TIEGCMHOME=Path/To/TIEGCM" to set environment variable.{RESET}\nEnter TIEGCM model directory: ')
    TIEGCMHOME = os.environ["TIEGCMHOME"]
# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(TIEGCMHOME, "tiegcmrun")

# Path to template .inp file.
INP_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.inp")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")

BENCHAMRKS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, 'benchmarks.json')


def get_mtime(file_path):
    """
    Get the 'mtime' data from the given file path.

    Parameters:
    file_path (str): The path to the file.

    Returns:
    list: A list containing the 'mtime' data, padded with zeros if necessary.
    """
    ds = xr.open_dataset(file_path)
    if 'mtime' in ds.variables:
        mtime_data = ds['mtime'].values
    mtime_arr = pad(mtime_data, [(0, 0), (0, max(4 - mtime_data.shape[1], 0))], mode='constant').tolist()
    return mtime_arr

def interp2d(variable, inlat, inlon, outlat, outlon):
    """
    Interpolates a 2D variable from input latitude and longitude grid to output latitude and longitude grid.

    Parameters:
    variable (ndarray): 2D array of the variable to be interpolated.
    inlat (ndarray): 1D array of input latitudes.
    inlon (ndarray): 1D array of input longitudes.
    outlat (ndarray): 1D array of output latitudes.
    outlon (ndarray): 1D array of output longitudes.

    Returns:
    ndarray: 2D array of the interpolated variable on the output grid.
    """

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
    """
    Interpolates a 3-dimensional variable from one set of levels and grid points to another set of levels and grid points.

    Args:
        variable (ndarray): The input variable with shape (ninlev, inlat, inlon).
        inlev (ndarray): The input levels.
        inlat (ndarray): The input latitudes.
        inlon (ndarray): The input longitudes.
        outlev (ndarray): The output levels.
        outlat (ndarray): The output latitudes.
        outlon (ndarray): The output longitudes.
        extrap (str): The extrapolation method. Must be one of 'constant', 'linear', or 'exponential'.

    Returns:
        ndarray: The interpolated variable with shape (noutlev, noutlat, noutlon).
    """
    
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
    """
    Interpolate and create a new primary file from an old TIEGCM primary file.

    Parameters:
    - fin (str): The filename of the old TIEGCM primary file.
    - hres (float): The horizontal resolution for the new primary file.
    - vres (float): The vertical resolution for the new primary file.
    - zitop (float): The top altitude for the new primary file.
    - fout (str): The filename of the new primary file to be created.

    Returns:
    None
    """
    print(f"Interpolating primary file {fin} to create new primary file {fout} at horizontal resolution {hres} and vertical resolution {vres} with zitop {zitop}.")
    # Some additional attributes for 4D fields
    lower_cap = 1e-12
    fill_top = ['TN', 'UN', 'VN', 'OP', 'TI', 'TE', 'N2D', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'OP_NM']
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
            'NO_NM': 'exponential', 'AR_NM': 'exponential'}

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

    # Longitude wrap is handled in interp2d, skip

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

        # The following variables never appear in standard TIEGCM runs
        # But they showed up in some non-standard TIEGCM primary histories
        # If that happens, skip those variables (the list may expand)
        elif varname in ['lat_bnds', 'lon_bnds', 'gw', 'area']:
            continue

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
            # In case HE is not in the old file (non-standard format), it has to be added to the new file
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



def segment_time(start_time_str, stop_time_str, interval_array):
    """
    Generate a list of time intervals between a start time and stop time based on a given interval array.

    Args:
        start_time_str (str): The start time in the format '%Y-%m-%dT%H:%M:%S'.
        stop_time_str (str): The stop time in the format '%Y-%m-%dT%H:%M:%S'.
        interval_array (list): A list containing the interval values in the order [days, hours, minutes, seconds].

    Returns:
        list: A list of time intervals in the format [[start_time, end_time], [start_time, end_time], ...].
    """
    # Convert start_time and stop_time to datetime objects
    start = datetime.strptime(start_time_str, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time_str, '%Y-%m-%dT%H:%M:%S')
    
    # Extract interval values from the array
    days, hours, minutes, seconds = interval_array
    
    # Create a timedelta object using the interval values
    delta = timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)
    
    # Generate the intervals
    intervals = []
    current = start
    while current < stop:
        next_time = min(current + delta, stop)
        intervals.append([current.strftime('%Y-%m-%dT%H:%M:%S'), next_time.strftime('%Y-%m-%dT%H:%M:%S')])
        current = next_time
    
    return intervals



def inp_pri_date(start_date_str, stop_date_str):
    """
    Convert start and stop date strings to datetime objects and extract relevant information.

    Args:
        start_date_str (str): Start date string in the format "%Y-%m-%dT%H:%M:%S".
        stop_date_str (str): Stop date string in the format "%Y-%m-%dT%H:%M:%S".

    Returns:
        tuple: A tuple containing the following values:
            - START_YEAR (int): The year of the start date.
            - START_DAY (int): The day of the year of the start date.
            - PRISTART (list): A list containing the day of the year, hour, minute, and second of the start date.
            - PRISTOP (list): A list containing the day of the year, hour, minute, and second of the stop date.
    """
    # Parse the start and stop dates
    start_time = datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
    stop_time = datetime.strptime(stop_date_str, "%Y-%m-%dT%H:%M:%S")

    # Extract START_YEAR and START_DAY
    START_YEAR = start_time.year
    START_DAY = start_time.timetuple().tm_yday

    # Format PRISTART and PRISTOP
    PRISTART = [start_time.timetuple().tm_yday, start_time.hour, start_time.minute, start_time.second]
    PRISTOP = [stop_time.timetuple().tm_yday, stop_time.hour, stop_time.minute, stop_time.second]

    return START_YEAR, START_DAY, PRISTART, PRISTOP

def valid_hist(start_time, stop_time):
    """
    Calculate valid divisions for a given time range.

    Parameters:
    start_time (str): The start date in the format '%Y-%m-%dT%H:%M:%S'.
    stop_time (str): The stop date in the format '%Y-%m-%dT%H:%M:%S'.

    Returns:
    list: A list of valid divisions for days, hours, minutes, and seconds.
    Each division is represented as a list [days, hours, minutes, seconds].
    """

    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    total_duration = stop - start
    total_seconds = total_duration.total_seconds()
    
    valid_divisions = []
    
    # Calculate valid divisions for days
    total_days = total_duration.days
    for n_day in range(1, total_days + 1):
        valid_divisions.append([n_day,0,0,0])
    
    # Calculate valid divisions for hours
    hours_divisions = [1, 2, 3, 4, 6, 12, 18, 24]
    for n_hour in hours_divisions:
        if total_seconds % (n_hour * 3600) == 0:
            valid_divisions.append([0,n_hour,0,0])
            
    # Calculate valid divisions for minutes
    minutes_divisions = [1, 2, 5, 10, 15, 30, 45, 60]
    for n_min in minutes_divisions:
        if total_seconds % (n_min * 60) == 0:
            valid_divisions.append([0,0,n_min,0])
    
    # Calculate valid divisions for seconds
    seconds_divisions = [1, 2, 5, 10, 15, 30, 45, 60]
    for n_sec in seconds_divisions:
        if total_seconds % n_sec == 0:
            valid_divisions.append([0,0,0,n_sec])
    
    return valid_divisions

def inp_mxhist(start_time, stop_time, x_hist, mxhist_warn):
    """
    Calculate the MXHIST value based on the given start and stop times, and x_hist values.

    Args:
        start_time (str): The start time in the format '%Y-%m-%dT%H:%M:%S'.
        stop_time (str): The stop time in the format '%Y-%m-%dT%H:%M:%S'.
        x_hist (tuple): A tuple containing the number of days, hours, minutes, and seconds for x_hist.
        mxhist_warn (str): A warning message for MXHIST.

    Returns:
        tuple: A tuple containing the calculated MXHIST value and the updated mxhist_warn message.

    Raises:
        None

    """
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    total_duration = stop - start
    total_seconds = total_duration.total_seconds()
    
    seconds_in_day = 86400
    seconds_in_hour = 3600
    seconds_in_min = 60
    
    n_day, n_hour, n_min, n_sec = x_hist
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec
    
    if step_seconds == 0:
        return "Invalid prihist: step cannot be 0."
    
    mxhist_day = seconds_in_day / step_seconds
    mxhist_hour = seconds_in_hour / step_seconds
    mxhist_min = seconds_in_min / step_seconds
    if mxhist_day >= 1:
        mxhist_warn = (mxhist_warn + "\n" if mxhist_warn is not None else "") + f"For a Daily output set MXHIST to {int(mxhist_day)}"
    if mxhist_hour >= 1:
        mxhist_warn = mxhist_warn +  f"\nFor a Hourly output set MXHIST to {int(mxhist_hour)}"
    if mxhist_min >= 1:
        mxhist_warn = mxhist_warn + f"\nFor a Minutely output set MXHIST to {int(mxhist_min)}"
    MXHIST = mxhist_day
    return(int(MXHIST), mxhist_warn)


def inp_sechist(SECSTART, SECSTOP):
    """
    Determines the value of SECHIST based on the given SECSTART and SECSTOP.

    Parameters:
    SECSTART (list): A list containing the start day.
    SECSTOP (list): A list containing the stop day.

    Returns:
    list: A list containing the value of SECHIST.

    """
    PRISTART_DAY = SECSTART[0]
    PRISTOP_DAY = SECSTOP[0]
    n_split_day = int(PRISTOP_DAY - PRISTART_DAY)
    if n_split_day >= 7:
        SECHIST = [1, 0, 0, 0]
    else:
        SECHIST = [0, 1, 0, 0]
    return SECHIST

def inp_prihist(PRISTART, PRISTOP):
    """
    Calculate the PRIHIST list based on the PRISTART and PRISTOP values.

    Parameters:
    PRISTART (list): A list containing the start day of the PRIHIST period.
    PRISTOP (list): A list containing the stop day of the PRIHIST period.

    Returns:
    list: The PRIHIST list, which is either [1, 0, 0, 0] or [0, 1, 0, 0] based on the number of days in the PRIHIST period.
    """
    PRISTART_DAY = PRISTART[0]
    PRISTOP_DAY = PRISTOP[0]
    n_split_day = int(PRISTOP_DAY - PRISTART_DAY)
    if n_split_day >= 7:
        PRIHIST = [1, 0, 0, 0]
    else:
        PRIHIST = [0, 1, 0, 0]
    return PRIHIST

def inp_pri_out(start_time, stop_time, PRIHIST, MXHIST_PRIM, pri_files, histdir, run_name):
    """
    Generate the output file names and the total number of files based on the given parameters.

    Parameters:
    start_time (str): The start time in the format '%Y-%m-%dT%H:%M:%S'.
    stop_time (str): The stop time in the format '%Y-%m-%dT%H:%M:%S'.
    PRIHIST (tuple): A tuple containing the number of days, hours, minutes, and seconds for PRIHIST.
    MXHIST_PRIM (int): The maximum number of primary history files.
    pri_files (int): The number of existing primary history files.
    histdir (str): The directory where the history files are stored.
    run_name (str): The name of the run.

    Returns:
    tuple: A tuple containing the output file names and the updated number of primary history files.

    """
    # Convert start and stop times to datetime
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    
    # Calculate total duration in seconds
    total_seconds = (stop - start).total_seconds()
    
    # Convert prihist to seconds
    n_day, n_hour, n_min, n_sec = PRIHIST
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec
    
    # Calculate model data per output file in seconds
    data_per_file_seconds = step_seconds * int(MXHIST_PRIM)
    
    # Calculate the total number of files, rounding up
    number_of_files = ceil(total_seconds / data_per_file_seconds)
    pri_files_n = pri_files + number_of_files
    if pri_files == 0:
        if number_of_files == 1:
                OUTPUT = OUTPUT = f"'{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc' , '{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files+1)}.nc'"
        else:
            PRIM_0 = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
            PRIM_N = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files_n)}.nc"
            OUTPUT = f"'{PRIM_0}','to','{PRIM_N}','by','1'"
    else:
        if number_of_files == 1:
            OUTPUT = OUTPUT = f"'{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc' , '{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files+1)}.nc'"
        else:
            PRIM_0 = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
            PRIM_N = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files_n)}.nc"
            OUTPUT = f"'{PRIM_0}','to','{PRIM_N}','by','1'"
    return OUTPUT, pri_files_n

def inp_sec_out(start_time, stop_time, SECHIST, MXHIST_SECH, sec_files, histdir, run_name):
    """
    Calculate the output file names and the total number of files for the secondary history output.

    Args:
        start_time (str): The start time of the simulation in the format '%Y-%m-%dT%H:%M:%S'.
        stop_time (str): The stop time of the simulation in the format '%Y-%m-%dT%H:%M:%S'.
        SECHIST (list): A list of integers representing the duration of each secondary history output file in days, hours, minutes, and seconds.
        MXHIST_SECH (int): The maximum number of secondary history output files per primary history output file.
        sec_files (int): The number of existing secondary history output files.
        histdir (str): The directory where the history output files are stored.
        run_name (str): The name of the simulation run.

    Returns:
        tuple: A tuple containing the output file name pattern and the updated number of secondary history output files.

    """
    # Convert start and stop times to datetime
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    sechist_delta = timedelta(days=SECHIST[0], hours=SECHIST[1], minutes=SECHIST[2], seconds=SECHIST[3])
    start = start + sechist_delta

    # Calculate total duration in seconds
    total_seconds = (stop - start).total_seconds()

    # Convert prihist to seconds
    n_day, n_hour, n_min, n_sec = SECHIST
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec

    # Calculate model data per output file in seconds
    data_per_file_seconds = step_seconds * int(MXHIST_SECH)

    # Calculate the total number of files, rounding up
    number_of_files = ceil(total_seconds / data_per_file_seconds)
    sec_files_n = sec_files + number_of_files

    if sec_files == 0:
        if number_of_files == 1:
            SECOUT = f"'{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files+1)}.nc'"
        else:
            SECH_0 = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files+1)}.nc"
            SECH_N = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_n+1)}.nc"
            SECOUT = f"'{SECH_0}','to','{SECH_N}','by','1'"
    else:
        if number_of_files == 1:
            SECOUT = f"'{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files+1)}.nc'"
        else:
            SECH_0 = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files+1)}.nc"
            SECH_N = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_n+1)}.nc"
            SECOUT = f"'{SECH_0}','to','{SECH_N}','by','1'"
    return SECOUT, sec_files_n

def inp_sec_date(start_time, stop_time, SECHIST):
    """
    Calculate the SECSTART and SECSTOP values based on the given start_time, stop_time, and SECHIST.

    Parameters:
    start_time (str): The start time in the format '%Y-%m-%dT%H:%M:%S'.
    stop_time (str): The stop time in the format '%Y-%m-%dT%H:%M:%S'.
    SECHIST (list): A list containing the number of days, hours, minutes, and seconds to be added to the start time.

    Returns:
    tuple: A tuple containing the SECSTART and SECSTOP values.

    Example:
    start_time = '2022-01-01T00:00:00'
    stop_time = '2022-01-02T00:00:00'
    SECHIST = [1, 0, 0, 0]
    inp_sec_date(start_time, stop_time, SECHIST)
    # Output: ([1, 0, 0, 0], [2, 0, 0, 0])
    """
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    sechist_delta = timedelta(days=SECHIST[0], hours=SECHIST[1], minutes=SECHIST[2], seconds=SECHIST[3])
    start = start + sechist_delta
    SECSTART = [start.timetuple().tm_yday,start.hour,start.minute,start.second]
    SECSTOP = [stop.timetuple().tm_yday,stop.hour,stop.minute,stop.second]

    return SECSTART, SECSTOP

def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--clobber", action="store_false",
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--onlycompile","-oc", action="store_true",
        help="Only Compile Tiegcm (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default=None,
        help="User mode (BASIC|INTERMEDIATE|EXPERT) (default: %(default)s)."
    )
    parser.add_argument(
        "--coupling","-co", action="store_true",
        help="Enable coupling (default: %(default)s)."
    )
    parser.add_argument(
        "--options_path", "-o", default=None,
        help="Path to JSON file of options (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "--execute","-e", action="store_true",
        help="Execute TIEGCM (default: %(default)s)."
    )
    parser.add_argument(
        "--compile","-c", action="store_true",
        help="Compile TIEGCM (default: %(default)s)."
    )
    parser.add_argument(
        "--benchmark", "-bench", default = None, type=valid_bench,
        help="Bechmark run name  (default: %(default)s)."
    )
    parser.add_argument(
        "--engage", default=None,
        help="Path to JSON file of engage options (default: %(default)s)."
    )
    return parser

def valid_bench(value):
    """
    Validate the benchmark option.

    Args:
        value (str): The benchmark option to validate.

    Returns:
        str: The validated benchmark option.

    Raises:
        argparse.ArgumentTypeError: If the value is not a valid benchmark option.
    """
    # Custom validation logic
    if value not in [None, 
                     'seasons', 'decsol_smax', 'decsol_smin', 'junsol_smax', 'junsol_smin','mareqx_smax', 'mareqx_smin', 'sepeqx_smax', 'sepeqx_smin',
                     'storms', 'dec2006_heelis_gpi', 'dec2006_weimer_imf', 'jul2000_heelis_gpi', 'jul2000_weimer_imf', 'nov2003_heelis_gpi', 'nov2003_weimer_imf', 'whi2008_heelis_gpi', 'whi2008_weimer_imf',
                     'climatology', 'climatology_smax', 'climatology_smin' 
                    ]:
        raise argparse.ArgumentTypeError(f"{value} is not a valid benchmark option.")
    return value

def get_run_option(name, description, mode="BASIC", skip_parameters=[]):
    """Prompt the user for a single run option.

    Prompt the user for a single run option. If no user input is provided,
    the default value is returned for the option. If valids are provided, the
    new value is compared against the valids, and rejected if not in the
    valids list.

    Parameters
    ----------
    name : str, default None
        Name of option
    description : dict, default None
        Dictionary of metadata for the option.
    mode : str
        User experience mode: "BASIC", "INTERMEDIATE", or "ADVANCED".

    Returns
    -------
    value_str : str
        Value of option as a string.

    Raises
    ------
    None
    """
    # Extract prompt, default, and valids.
    level = description["LEVEL"]
    prompt = description.get("prompt", "")
    default = description.get("default", None)
    valids = description.get("valids", None)
    var_description = description.get("description", None)
    warning = description.get("warning", None)
    # Compare the current mode to the parameter level setting. If the variable
    # level is higher than the user mode, just use the default.
    fourvar_variables = ["SOURCE_START","segment","PRISTART","PRISTOP","PRIHIST","SECSTART","SECSTOP","SECHIST"]

    if mode == "BENCH" and level in ["BASIC","INTERMEDIATE", "EXPERT"]:
        if name in fourvar_variables and default is not None:
            return  ' '.join(map(str, default))
        else:
            return default
    if mode == "BASIC" and name in skip_parameters:
        if name in fourvar_variables and default is not None:
            return  ' '.join(map(str, default))
        else:
            return default    
    if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
        if name in fourvar_variables:
            return  ' '.join(map(str, default))
        else:
            return default
    if mode == "INTERMEDIATE" and level in ["EXPERT"]:
        if name in fourvar_variables:
            return  ' '.join(map(str, default))
        else:
            return default

    if warning is not None:
        print(f'{YELLOW}{warning}{RESET}')
    og_prompt = prompt
    file_variables = ["SOURCE","GPI_NCFILE","IMF_NCFILE","AMIENH","AMIESH","GSWM_MI_DI_NCFILE","GSWM_MI_SDI_NCFILE","GSWM_NM_DI_NCFILE","GSWM_NM_SDI_NCFILE","HE_COEFS_NCFILE","BGRDDATA_NCFILE","CTMT_NCFILE","SABER_NCFILE","TIDI_NCFILE","MIXFILE"]
    valid_bool = False
    bool_True = ["YES","Yes","yes","Y","y","TRUE","True","true","T","t","1",True,1]
    bool_False = ["NO","No","no","N","n","FALSE","False","false","F","f","0",False,0]
    array_variables = ["TIDE","TIDE2","NUDGE_FLDS","NUDGE_SPONGE","NUDGE_DELTA","NUDGE_POWER"]
    # If provided, add the valid values in val1|val2 format to the prompt.
    if valids is not None: 
        if name == "vertres":
            vs = "|".join(map(lambda x: str(Fraction(x)), valids))
            prompt += f" ({vs})"            
        else:
            if True in valids and False in valids:
                valid_bool = True
            vs = "|".join(map(str,valids))
            prompt += f" ({vs})"

    # If provided, add the default to the prompt.
    if default is not None:
        if name == "vertres":
            prompt += f" [{GREEN}{str(Fraction(default))}{RESET}]"
        else:
            prompt += f" [{GREEN}{default}{RESET}]"

    # Prompt the user and fetch the input until a good value is provided.
    ok = False
    option_value = ""
    first_pass_modules = True
    while not ok:
        if name in ("other_input", "other_pbs","other_job","job_chain"):
            prompt = og_prompt
            prompt += f" [{GREEN}{default}{RESET}]"
            temp_value = input(f"{prompt} / ENTER to go next: ")
            if temp_value != "":
                default.append(temp_value)
                option_value = default
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                option_value = default
                ok = True
            elif temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
        elif name == "SECFLDS":
            prompt = og_prompt
            prompt += f" [{GREEN}{default}{RESET}]"
            temp_value = input(f"{prompt} / ENTER to go next: ")
            if temp_value != "":
                if "," in temp_value:
                    temp_value = temp_value.replace("'", "")
                    default.extend(s.replace(" ", "") for s in temp_value.split(',')) 
                else:
                    temp_value = temp_value.replace("'", "")
                    default.extend(s.replace(" ", "") for s in temp_value.split()) 
                option_value = default
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                option_value = default
                ok = True
            elif temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
            """
            elif name in array_variables:
                prompt = og_prompt
                prompt += f" [{GREEN}{default}{RESET}]"
                temp_value = input(f"{prompt} / ENTER to go next: ")
                if temp_value != "":
                    if "," in temp_value:
                        temp_value = temp_value.replace("'", "")
                        default.extend(s.replace(" ", "") for s in temp_value.split(',')) 
                    else:
                        temp_value = temp_value.replace("'", "")
                        default.extend(s.replace(" ", "") for s in temp_value.split()) 
                    option_value = default
                elif temp_value == 'none' or temp_value == 'None':
                    option_value = json.loads('[null]')
                    ok = True
                elif temp_value == "":
                    option_value = default
                    ok = True
                elif temp_value == "?":
                    print(f'{YELLOW}{var_description}{RESET}')
            """
        elif name == "modules":
            prompt = og_prompt
            prompt += f" [{GREEN}{default}{RESET}]"
            temp_value = input(f"{prompt} / ENTER to go next: ")
            if temp_value != "":
                if first_pass_modules == True:
                    default = []
                    first_pass_modules = False
                if "\n" in temp_value:
                    temp_value = temp_value.replace("'", "")
                    temp_array = temp_value.split('\n')
                    temp_array = [s.replace("module load", "") for s in temp_array]
                    temp_array = [s.replace(" ", "") for s in temp_array]
                    default.extend(temp_array) 
                if "," in temp_value:
                    temp_value = temp_value.replace("'", "")
                    temp_array = temp_value.split(',')
                    temp_array = [s.replace("module load", "") for s in temp_array]
                    temp_array = [s.replace(" ", "") for s in temp_array]
                    default.extend(temp_array) 
                else:
                    temp_value = temp_value.replace("'", "")
                    temp_array = temp_value.split()
                    temp_array = [s.replace("module load", "") for s in temp_array]
                    temp_array = [s.replace(" ", "") for s in temp_array]
                    default.extend(temp_array)
                modules_array_cleaned = []
                for module in default:
                    if module not in ['module', 'load']:
                        modules_array_cleaned.append(module) 
                default = modules_array_cleaned
                option_value = default
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                option_value = default
                ok = True
            elif temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
        elif name in file_variables:
            option_value = input(f"{prompt}: ")
            # Use the default if no user input provided.
            if option_value == "":
                option_value = default
                if name == "SOURCE" and option_value == None:
                    continue
                else:
                    ok = True
            # Use None if user input is none or None.
            elif option_value == 'none' or option_value == 'None':
                option_value = None
                ok = True
            elif option_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
                continue
            elif option_value != None:
                if os.path.isfile(option_value) == False:
                    if os.path.isdir(option_value) == True:
                        print(f'{YELLOW} {option_value} is a directory. Please provide a file path.{RESET}')
                        continue
                    else:
                        file_path = find_file(option_value,TIEGCMDATA)
                        if file_path == None:
                            print(f'{YELLOW} Unable to find {option_value} in {TIEGCMDATA}.\n Give path to file as an alternative.{RESET}')
                            continue
                        else:
                            print(f'File Found: {file_path}')
                            option_value = str(file_path)
                            ok = True
                else:
                    option_value = str(option_value)
                    ok = True
        elif name in fourvar_variables:
            prompt = og_prompt
            if valids is not None: 
                prompt += " Example:"
                vs =' | '.join([','.join(map(str, sublist)) for sublist in valids])
                prompt += f" ({vs})"
            if default not in [None, [None]] :    
                default_print = ' '.join(map(str, default))
                prompt += f" [{GREEN}{default_print}{RESET}]"
            temp_value = input(f"{prompt}: ")
            temp_array = []
            if temp_value not in ["","?","none","None"]:
                if "," in temp_value:
                    temp_value = temp_value.replace("'", "")
                    temp_array.extend(s.replace(" ", "") for s in temp_value.split(',')) 
                    temp_array = [int(i) for i in temp_array]
                else:
                    temp_value = temp_value.replace("'", "")
                    temp_array.extend(s.replace(" ", "") for s in temp_value.split()) 
                    temp_array = [int(i) for i in temp_array]
                option_value = temp_array
                if len(option_value) != 4 or option_value == [0,0,0,0]:
                    print(f'{YELLOW}Invalid Value: {option_value}{RESET}')
                    continue
                else:
                    if valids is not None: 
                        if option_value not in valids:
                            print(f'{YELLOW}{option_value} not in default list. \nSetting dependent defaults/suggested values to None.{RESET}')
                            option_value =' '.join(map(str, option_value))
                            ok = True
                        else:
                            option_value =' '.join(map(str, option_value))
                            ok = True
                    else:
                        option_value =' '.join(map(str, option_value))
                        ok = True
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                if default not in [None,[None]]:
                    option_value = ' '.join(map(str, default))
                else:
                    option_value = [None]
                ok = True
            elif temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
        else:
            # Fetch input from the user.
            option_value = input(f"{prompt}: ")

            # Use the default if no user input provided.
            if option_value == "":
                option_value = default
            # Use None if user input is none or None.
            elif option_value == 'none' or option_value == 'None':
                option_value = None
            elif option_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
                continue
            # Validate the result. If bad, start over.
            if name == "vertres":
                if valids is not None and float(Fraction(option_value)) not in valids:
                    print(f"Invalid value for option {name}: {option_value}!")
                    continue
                else:
                    option_value=float(Fraction(option_value))
            elif name in ["horires","ELECTRON_HEATING","zitop","mres"]:
                if valids is not None and float(option_value) not in valids:
                    print(f"Invalid value for option {name}: {option_value}!")
                    continue
            elif name in ["AURORA","DYNAMO","CALC_HELIUM","EDDY_DIF","JOULEFAC","COLFAC","OPDIFFCAP","CURRENT_PG","CURRENT_KQ","NUDGE_ALPHA"]:
                if option_value != None: 
                    if float(option_value) not in valids:
                        print(f"Invalid value for option {name}: {option_value}!")
                        continue
            else:
                if valid_bool == True:
                    if option_value not in bool_True and option_value not in bool_False:
                        if mode !="EXPERT":
                            print(f"Invalid value for option {name}: {option_value}!")
                            continue
                elif valids is not None and option_value not in valids:
                    print(f"Invalid value for option {name}: {option_value}!")
                    continue

            # Keep this result.
            ok = True
            if valid_bool == True:
                if option_value in bool_True:
                    option_value = True
                elif option_value in bool_False:
                    option_value = False
            if option_value != None and type(option_value) != bool:
                option_value = str(option_value)
    # Return the option as a string.
    return option_value

def select_source_defaults(options, option_descriptions):
    """
    Select the default values for the 'source' option based on the given input options.

    Args:
        options (dict): A dictionary containing the input options.
        option_descriptions (dict): A dictionary containing the descriptions of the available options.

    Returns:
        str: The default value for the 'source' option.

    """
    start_time = options["inp"]["start_time"]
    time_dhms = time_to_dhms(start_time)
    flux_level = options["inp"]["solar_flux_level"]
    if flux_level == "low":
        f107 = 70
    elif flux_level == "medium":
        f107 = 140
    elif flux_level == "high":
        f107 = 200
    if time_dhms[0] >= 1 and time_dhms[0] <81:
        source_default = find_file(f"decsol_f{f107}*",TIEGCMDATA)
    elif time_dhms[0] >= 81 and time_dhms[0] <173:
        source_default = find_file(f'mareqx_f{f107}*',TIEGCMDATA)
    elif time_dhms[0] >= 173 and time_dhms[0] <265:
        source_default = find_file(f'junsol_f{f107}*',TIEGCMDATA)
    elif time_dhms[0] >= 265 and time_dhms[0] <356:
        source_default = find_file(f'seqex_f{f107}*',TIEGCMDATA)
    elif time_dhms[0] >= 356:
        source_default = find_file(f'decsol_f{f107}*',TIEGCMDATA)       
    return source_default    

def select_resource_defaults(options, option_descriptions):
    """
    Selects the default values for the 'select', 'ncpus', and 'mpiprocs' options based on the given input options.

    Args:
        options (dict): A dictionary containing the input options.
        option_descriptions (dict): A dictionary containing the descriptions of the available options.

    Returns:
        tuple: A tuple containing the default values for 'select', 'ncpus', and 'mpiprocs' options.

    """
    horires = options["model"]["specification"]["horires"]
    hpc_platform = options["simulation"]["hpc_system"]
    od = option_descriptions["job"][hpc_platform]
    o = options["job"]
    if hpc_platform == "derecho":
        od=od["resource"]
        for on in od:
            if on == "select":
                if float(horires) == 2.5 or float(horires) == 5:
                    select_default = 3
                elif float(horires) == 1.25:
                    select_default = 3
                elif float(horires) == 0.625:
                    select_default = 3         
            if on == "ncpus":
                if float(horires) == 2.5 or float(horires) == 5:
                    ncpus_default = 128
                elif float(horires) == 1.25:
                    ncpus_default = 128
                elif float(horires) == 0.625:
                    ncpus_default = 128
            if on == "mpiprocs":
                if float(horires) == 2.5 or float(horires) == 5:
                    mpiprocs_default = 96
                elif float(horires) == 1.25:
                    mpiprocs_default = 96
                elif float(horires) == 0.625:
                    mpiprocs_default = 96
    elif hpc_platform == "pleiades":
        od=od["resource"]
        o=o["resource"]
        if o["model"] == "bro":
            max_ncpus = 24
        elif o["model"] == "has":
            max_ncpus = 24
        elif o["model"] == "ivy":
            max_ncpus = 18
        elif o["model"] == "san":
            max_ncpus = 12
        for on in od:
            if on == "select":
                if float(horires) == 2.5 or float(horires) == 5:
                    select = 72/max_ncpus
                if float(horires) == 1.25:
                    select = 144/max_ncpus
                if float(horires) == 0.625:
                    select = 288/max_ncpus
                select_default = int(select)
            if on == "ncpus":
                ncpus_default = max_ncpus
            if on == "mpiprocs":
                mpiprocs_default = ncpus_default
    return select_default,ncpus_default,mpiprocs_default


def prompt_user_for_run_options(args):
    """Prompt the user for run options.

    Prompt the user for run options.

    NOTE: In this function, the complete set of parameters is split up
    into logical groups. This is done partly to make the organization of the
    parameters more obvious, and partly to allow the values of options to
    depend upon previously-specified options.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options

    Returns
    -------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Raises
    ------
    None
    """

    global TIEGCMDATA
    global TIEGCMHOME
    input_build_skip = False
    pbs_build_skip = False
    base_skip = False
    skip_parameters = []
    # Save the user mode.
    mode = args.mode
    benchmark = args.benchmark
    engage= args.engage
    coupling = args.coupling
    if engage != None:
        skip_parameters = engage["skip"]
    if benchmark != None and mode == None:
        mode = "BENCH"
    elif mode == None:
        mode = "BASIC"
    onlycompile = args.onlycompile
    if onlycompile == True:
        input_build_skip = True
        pbs_build_skip = True
        skip_parameters = ['input_file','log_file','job_name','modeldir','parentdir','tgcmdata']
        if coupling == True:
            skip_parameters.append('modelexe')
    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)
    with open(BENCHAMRKS_FILE, "r", encoding="utf-8") as f:
        benchmarks_options = json.load(f)
    # Initialize the dictionary of program options.
    options = {}

    #-------------------------------------------------------------------------
    if base_skip == False:
        # General options for the simulation
        o = options["simulation"] = {}
        od = option_descriptions["simulation"]
        if benchmark != None:
            od["job_name"]["default"] = benchmark
        elif engage != None:
            od["job_name"]["default"] = engage["job_name"]
        system_name = os.popen('hostname').read().strip()
        if 'pfe' in system_name.lower():
            od["hpc_system"]["default"]= "pleiades"
        elif 'derecho' in system_name.lower():
            od["hpc_system"]["default"] = "derecho"
        else:
            od["hpc_system"]["default"] = "linux"
        # Prompt for the parameters.
        for on in ["job_name", "hpc_system"]:
            o[on] = get_run_option(on, od[on], mode, skip_parameters)
    #-------------------------------------------------------------------------

    # Model options
    options["model"] = {}
    #------------------------------------

    # Data Options
     
    options["model"]["data"] = {}
    o = options["model"]["data"]
    temp_mode = mode
    od = option_descriptions["model"]["data"]
    od["modeldir"]["default"] = TIEGCMHOME 
    o["modeldir"] = get_run_option("modeldir", od["modeldir"], mode, skip_parameters)
    if engage != None:
        od["parentdir"]["default"] = engage["parentdir"]
        o["parentdir"] = get_run_option("parentdir", od["parentdir"], mode, skip_parameters)
        od["execdir"]["default"] = o["parentdir"]
        od["workdir"]["default"] = o["parentdir"]
        od["histdir"]["default"] = o["parentdir"]
        o["execdir"] = get_run_option("execdir", od["execdir"], mode, skip_parameters)
    else:
        od["parentdir"]["default"] = os.getcwd() 
        o["parentdir"] = get_run_option("parentdir", od["parentdir"], mode, skip_parameters)
        if o["parentdir"] == None:
            temp_mode = "INTERMEDIATE"
            o["execdir"] = get_run_option("execdir", od["execdir"], temp_mode, skip_parameters)
            od["workdir"]["default"] = o["execdir"]
            od["histdir"]["default"] = o["execdir"]
        else:
            od["execdir"]["default"] = os.path.join(o["parentdir"],"exec")
            od["workdir"]["default"] = os.path.join(o["parentdir"],"stdout")
            od["histdir"]["default"] = os.path.join(o["parentdir"],"hist")
            o["execdir"] = get_run_option("execdir", od["execdir"], mode, skip_parameters)

    o["workdir"] = get_run_option("workdir", od["workdir"], temp_mode, skip_parameters)
    o["histdir"] = get_run_option("histdir", od["histdir"], temp_mode, skip_parameters)
    temp_mode = mode    
    o["utildir"] = os.path.join(o["modeldir"],'scripts')
    od["tgcmdata"]["default"] = TIEGCMDATA
    o["tgcmdata"] = get_run_option("tgcmdata", od["tgcmdata"], mode, skip_parameters)
    if base_skip == False:
        if benchmark == None or engage == None:
            o["input_file"] = get_run_option("input_file", od["input_file"], mode, skip_parameters)
            if o["input_file"]  != None:
                input_build_skip = True
        else:
            o["input_file"] = None
    od["log_file"]["default"] =  f'{o["workdir"]}/{options["simulation"]["job_name"]}.out'
    o["log_file"] = get_run_option("log_file", od["log_file"], mode, skip_parameters)

    if od["make"]["default"] == None:
        if options["simulation"]["hpc_system"] == "derecho":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_de')
        elif options["simulation"]["hpc_system"] == "pleiades":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_pf')
        elif options["simulation"]["hpc_system"] == "linux":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_linux')
    o["make"] = get_run_option("make", od["make"], mode, skip_parameters)
    od["modelexe"]["default"] = os.path.join(o["execdir"],"tiegcm.exe")
    o["modelexe"] = get_run_option("modelexe", od["modelexe"], mode, skip_parameters)
    if os.path.isfile(o["modelexe"]) == False:
        o["modelexe"] = os.path.join(o["execdir"],o["modelexe"])
        if args.compile == False and args.onlycompile == False:
            print(f'{YELLOW}Unable to find {o["modelexe"]}, model must be compiled. Use --compile/-c or --onlycompile/-oc {RESET}')
    if args.coupling == True:
        od["coupled_modelexe"]["default"] = os.path.join(o["execdir"],"tiegcm.x")
        o["coupled_modelexe"] = get_run_option("coupled_modelexe", od["coupled_modelexe"], mode, skip_parameters)
        if os.path.isfile(o["coupled_modelexe"]) == False:
            o["coupled_modelexe"] = os.path.join(o["execdir"],o["coupled_modelexe"])
            #o["modelexe"] = o["coupled_modelexe"]
            if args.compile == False and args.onlycompile == False:
                print(f'{YELLOW}Unable to find {o["coupled_modelexe"]}, model must be compiled. Use --compile/-c or --onlycompile/-oc {RESET}')

    TIEGCMDATA = o["tgcmdata"]
    
    #------------------------------------

    # Specification Options
            
    options["model"]["specification"] = {}
    o = options["model"]["specification"]

    od = option_descriptions["model"]["specification"]
                   
    for on in od:
        if engage != None:
            od["horires"]["default"] = engage["horires"]
        if on == "vertres":
            od["vertres"]["default"] = vertres
        elif on == "mres":
            od["mres"]["default"] = mres
        o[on] = get_run_option(on, od[on], mode, skip_parameters)
        if on =="horires":
            horires = float(o[on])
            vertres, mres, nres_grid, STEP = resolution_solver(horires)

    if o["nres_grid"] == None:
        o["nres_grid"] = nres_grid

    
    if input_build_skip == True:
        input_build_skip = True
        o["segmentation"] = False
    #-------------------------------------------------------------------------
    # INP options
    if input_build_skip == False:        

        options["inp"] = {}
        o = options["inp"]

        od = option_descriptions["inp"]     
        od["STEP"]["default"] = STEP

        run_name = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
        histdir = options["model"]["data"]["histdir"]
        if benchmark != None:
            oben = benchmarks_options[benchmark]["inp"]
            for on in od:
                if on in oben:
                    if on == "SOURCE":
                        od[on]["default"] = find_file('*tiegcm*'+options['simulation']['job_name']+'*.nc', TIEGCMDATA)
                    elif on in ["OUTPUT", "SECOUT"]:
                        temp_output = oben[on]
                        try:
                            temp_output = temp_output.replace("+histdir+", histdir)
                            temp_output = temp_output.replace("+run_name+", run_name)
                        except:
                            temp_output = temp_output
                        od[on]["default"] = temp_output
                    elif on in ["GPI_NCFILE","IMF_NCFILE"]:
                        od[on]["default"] = find_file(oben[on], TIEGCMDATA)
                    elif on in ["other_input"]:
                        temp_output = [item.replace("+tiegcmdata+", TIEGCMDATA) if item is not None and item != 'null' else item for item in oben[on]]
                        od[on]["default"] = temp_output
                    elif oben[on] == None:
                        od[on]["default"] = od[on]["default"]
                    else:
                        od[on]["default"] = oben[on]
                
        od["LABEL"]["default"] = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
        if options["simulation"]["hpc_system"] == "pleiades":
            od["SECFLDS"]["warning"] = "Limit SECFLDS. File libraries on pleiades are built without big file support."
        temp_mode = mode
        skip_inp = []
        start_stop_set = 0
        for on in od:
            if start_stop_set == 0 and benchmark == None:
                temp_mode =  "BASIC"  # "INTERMEDIATE"
            if on == "start_time" and benchmark != None:
                continue
            elif on == "stop_time" and benchmark != None:
                continue
            elif on == "start_time" and benchmark == None:
                if engage != None:
                    od["start_time"]["default"] = engage["start_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "stop_time" and benchmark == None:
                if engage != None:
                    od["stop_time"]["default"] = engage["stop_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] != None:
                    start_stop_set = 1
                    temp_mode = mode
                    START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(o["start_time"], o["stop_time"])
                    od["START_YEAR"]["default"] = START_YEAR
                    od["START_DAY"]["default"] = START_DAY
                    od["PRISTART"]["default"] = PRISTART
                    od["PRISTOP"]["default"] = PRISTOP
            elif on == "secondary_start_time" and benchmark == None:
                od["secondary_start_time"]["default"] = o["start_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                #od["SECSTART"]["default"] = o[on]
            elif on == "secondary_stop_time" and benchmark == None:
                od["secondary_stop_time"]["default"] = o["stop_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                #od["SECSTOP"]["default"] = o[on]
            elif on == "segment" and benchmark == None:
                if engage != None:
                    od["segment"]["default"] = engage["segment"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] != [None]:
                    options["model"]["specification"]["segmentation"] = True
                    segment = [int(i) for i in o[on].split()]
                    runtimes = segment_time(o["start_time"], o["stop_time"], segment)
                    segment_warn_0 = f"Segmentation is set to {segment}.\n"
                    segment_warn_1 = f" is set for one segment\neg.{runtimes[0][0]} to {runtimes[0][1]}" 
                    od["PRIHIST"]["warning"] = (od["PRIHIST"]["warning"] + "\n" if od["PRIHIST"]["warning"]  is not None else "") + segment_warn_0 + "PRIHIST" + segment_warn_1
                    od["MXHIST_PRIM"]["warning"] = (od["MXHIST_PRIM"]["warning"] + "\n" if od["MXHIST_PRIM"]["warning"]  is not None else "") + segment_warn_0 + "MXHIST_PRIM" + segment_warn_1
                    od["SECHIST"]["warning"] = (od["SECHIST"]["warning"] + "\n" if od["SECHIST"]["warning"]  is not None else "") + segment_warn_0 + "SECHIST" + segment_warn_1
                    od["MXHIST_SECH"]["warning"] = (od["MXHIST_SECH"]["warning"] + "\n" if od["MXHIST_SECH"]["warning"]  is not None else "") + segment_warn_0 + "MXHIST_SECH" + segment_warn_1
                    od["OUTPUT"]["warning"] = (od["OUTPUT"]["warning"] + "\n" if od["OUTPUT"]["warning"]  is not None else "") + "Primary Output can be ignored. Will be set on segmentation"
                    od["SECOUT"]["warning"] = (od["SECOUT"]["warning"] + "\n" if od["SECOUT"]["warning"]  is not None else "") + "Secondary Output can be ignored. Will be set on segmentation"
            elif on == "solar_flux_level" and benchmark == None:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "SOURCE":
                if benchmark == None:
                    od["SOURCE"]["default"] = select_source_defaults(options, option_descriptions)    
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] == None:
                    o[on] = get_run_option(on, od[on], "BASIC")
            elif on == "SOURCE_START":
                od["SOURCE_START"]["valids"] = get_mtime(options["inp"]["SOURCE"])
                #for arr in mtime_arr:
                temp_mode_1 = temp_mode
                if len(od["SOURCE_START"]["valids"]) > 1:    
                    temp_mode_1 = "INTERMEDIATE"
                od["SOURCE_START"]["default"] = od["SOURCE_START"]["valids"][0]
                o[on] = get_run_option(on, od[on], temp_mode_1)
            elif on == "PRIHIST" and benchmark== None:
                PRIHIST = inp_prihist(PRISTART,PRISTOP)
                od["PRIHIST"]["default"] =  PRIHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                PRIHIST = [int(i) for i in o[on].split()]
                MXHIST_PRIM_set ,MXHIST_PRIM_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], PRIHIST, od["MXHIST_PRIM"]["warning"])
                od["MXHIST_PRIM"]["default"] = MXHIST_PRIM_set
                od["MXHIST_PRIM"]["warning"] = MXHIST_PRIM_warning_set
            elif on == "MXHIST_PRIM" and benchmark== None:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                MXHIST_PRIM = int(o[on])
            elif on == "OUTPUT" and benchmark== None:
                OUTPUT, pri_files_n = inp_pri_out(o["start_time"], o["stop_time"], PRIHIST, MXHIST_PRIM, 0, histdir,run_name)
                od["OUTPUT"]["default"] = OUTPUT
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "SECHIST" and benchmark== None:
                SECHIST = inp_sechist(PRISTART,PRISTOP)
                od["SECHIST"]["default"] =  SECHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                SECHIST = [int(i) for i in o[on].split()]
                MXHIST_SECH_set ,MXHIST_SECH_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], SECHIST, od["MXHIST_SECH"]["warning"])
                od["MXHIST_SECH"]["default"] = MXHIST_SECH_set
                od["MXHIST_SECH"]["warning"] = MXHIST_SECH_warning_set
                SECSTART, SECSTOP = inp_sec_date(o["secondary_start_time"], o["secondary_stop_time"], SECHIST)
                od["SECSTART"]["default"] = SECSTART
                od["SECSTOP"]["default"] = SECSTOP                
            elif on == "MXHIST_SECH" and benchmark== None:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                MXHIST_SECH = int(o[on])
            elif on == "SECOUT" and benchmark== None:
                SECOUT, sec_files_n = inp_sec_out(o["secondary_start_time"], o["secondary_stop_time"],  SECHIST, MXHIST_SECH, 0, histdir,run_name)
                od["SECOUT"]["default"] = SECOUT
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "POTENTIAL_MODEL":
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] == "HEELIS":
                    skip_inp_temp = ["IMF_NCFILE","BXIMF","BYIMF","BZIMF","SWDEN","SWVEL"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
                elif o[on] == "WEIMER":
                    skip_inp_temp = ["CTPOTEN"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on == "ONEWAY":
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)            
            elif on == "GPI_NCFILE" and on not in skip_inp:
                if engage != None:
                    od["GPI_NCFILE"]["default"] =  find_file('gpi_*', TIEGCMDATA)
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] != None:
                    skip_inp_temp = ["KP","POWER","CTPOTEN","F107","F107A"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
                    od["F107"]["warning"] = "F10.7 can be read by GPI File and can be skipped."
                    od["F107A"]["warning"] = "81-Day Average of F10.7 can be read by GPI File and can be skipped."
            elif on == "IMF_NCFILE" and on not in skip_inp:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] != None:
                    skip_inp_temp = ["BXIMF","BYIMF","BZIMF","SWDEN","SWVEL"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on == "KP" and on not in skip_inp:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                if o[on] != None:
                    skip_inp_temp = ["POWER","CTPOTEN"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on == "GSWM_MI_DI_NCFILE":
                GSWM_MI_DI_NCFILE = find_file(f'*gswm_diurn_{horires}d_99km*', TIEGCMDATA)
                od[on]["default"] = f"{GSWM_MI_DI_NCFILE}" if GSWM_MI_DI_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "GSWM_MI_SDI_NCFILE":
                GSWM_MI_SDI_NCFILE = find_file(f'*gswm_semi_{horires}d_99km*', TIEGCMDATA)
                od[on]["default"] = f"{GSWM_MI_SDI_NCFILE}" if GSWM_MI_SDI_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "GSWM_NM_DI_NCFILE":
                GSWM_NM_DI_NCFILE = find_file(f'*gswm_nonmig_diurn_{horires}d_99km*', TIEGCMDATA)
                od[on]["default"] = f"{GSWM_NM_DI_NCFILE}" if GSWM_NM_DI_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "GSWM_NM_SDI_NCFILE":
                GSWM_NM_SDI_NCFILE = find_file(f'*gswm_nonmig_semi_{horires}d_99km*', TIEGCMDATA)
                od[on]["default"] = f"{GSWM_NM_SDI_NCFILE}" if GSWM_NM_SDI_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on == "HE_COEFS_NCFILE":
                HE_COEFS_NCFILE = f"{find_file(f'*he_coefs_dres*', TIEGCMDATA)}"
                od[on]["default"] = f"{HE_COEFS_NCFILE}" if HE_COEFS_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on not in skip_inp:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
            elif on in skip_inp:
                o[on] = od[on]["default"]
            

    #-------------------------------------------------------------------------
    hpc_platform = options["simulation"]["hpc_system"]
    if hpc_platform == "linux":
        pbs_build_skip = True
    if pbs_build_skip == False:
        # PBS options
        options["job"] = {}
        o = options["job"]
        skip_pbs = []
        # Common (HPC platform-independent) options
        od = option_descriptions["job"]["_common"]
        od["account_name"]["default"] = os.getlogin()
        for on in od:
            if engage != None:
                    od['account_name']['default'] = engage['account_name']
            o[on] = get_run_option(on, od[on], mode, skip_parameters)

        # HPC platform-specific options
        hpc_platform = options["simulation"]["hpc_system"]
        od = option_descriptions["job"][hpc_platform]
        if hpc_platform == "derecho":
            o["mpi_command"] = "mpirun"
            if engage != None:
                    od['queue']['default'] = engage['queue']
                    od['job_priority']['default'] = engage['job_priority']
                    od['walltime']['default'] = engage['walltime']
        elif hpc_platform == "pleiades":
            o["mpi_command"] = "mpiexec_mpt"
        for on in od:
            if on == "resource":
                options["job"]["resource"] = {}
                odt = od["resource"]
                ot = options["job"]["resource"]
                for ont in odt:
                    if hpc_platform == "derecho":
                        select_default,ncpus_default,mpiprocs_default = select_resource_defaults(options,option_descriptions)
                        odt["select"]["default"] = select_default
                        odt["ncpus"]["default"] = ncpus_default
                        odt["mpiprocs"]["default"] = mpiprocs_default
                    elif hpc_platform == "pleiades":
                        if ont == "model":
                            ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
                            select_default,ncpus_default,mpiprocs_default = select_resource_defaults(options,option_descriptions)
                            odt["select"]["default"] = select_default
                            odt["ncpus"]["default"] = ncpus_default
                            odt["mpiprocs"]["default"] = mpiprocs_default
                    if ont == "select":
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
                        nnodes = ot[ont]
                    elif ont == "ncpus":
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
                        ncpus = ot[ont]
                    elif ont == "mpiprocs":
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
                        mpiprocs = ot[ont]
                    else:
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
            elif on =="nprocs":
                od[on]["default"] = int(nnodes) * int(mpiprocs) #int(mpiprocs)
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
            elif on == "module_list":
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
                if o[on] != None:
                    skip_pbs.append("modules")
            elif on == "modules" and on not in skip_pbs:
                if engage != None:
                    od[on]["default"] = engage["modules"]
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
            elif on == "project_code":
                if engage != None:
                    od[on]["default"] = engage["project_code"]
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
            elif on not in skip_pbs:
                o[on] = get_run_option(on, od[on], mode, skip_parameters)

    # Return the options dictionary.
    return options

def resolution_solver(horires, engage_options=None):
    if float(horires) == 5:
        vertres = 0.5
        mres = 2
        STEP = 60
    elif float(horires) == 2.5:
        vertres = 0.25
        mres = 2
        STEP = 30
    elif float(horires) == 1.25:
        vertres = 0.125
        mres = 1
        STEP = 10
    elif float(horires) == 0.625:
        vertres = 0.0625
        mres = 0.5
        STEP = 5

    if mres == 2:
        nres_grid = 5
    elif mres == 1:
        nres_grid = 6
    elif mres == 0.5:
        nres_grid = 7
    
    if engage_options != None:
        STEP = engage_options["STEP"]
    
    return vertres, mres, nres_grid, STEP



def compile_tiegcm(options, debug, coupling):
    """
    Compiles the TIEGCM model with the given options.

    Args:
        options (dict): A dictionary containing the model options.
        debug (bool): A boolean indicating whether to enable debug mode.
        coupling (bool): A boolean indicating whether to enable coupling.

    Returns:None
    """
    o = options
    modeldir  = o["model"]["data"]["modeldir"]
    execdir   = o["model"]["data"]["execdir"]
    workdir = o["model"]["data"]["workdir"]
    outdir = o["model"]["data"]["histdir"]
    tgcmdata  = o["model"]["data"]["tgcmdata"]
    utildir   = os.path.join(o["model"]["data"]["modeldir"],"scripts")
    try:
        input     = o["model"]["data"]["input_file"]
    except:
        input = ""
    try:
        output    = o["model"]["data"]["log_file"]
    except:
        output = ""
    horires   = float(o["model"]["specification"]["horires"])
    vertres   = float(o["model"]["specification"]["vertres"])
    zitop     = float(o["model"]["specification"]["zitop"])
    mres      = float(o["model"]["specification"]["mres"])
    nres_grid = float(o["model"]["specification"]["nres_grid"])
    make      = o["model"]["data"]["make"]
    coupling  = coupling
    if coupling == True:
        modelexe = os.path.basename(o["model"]["data"]["coupled_modelexe"])
        model = o["model"]["data"]["coupled_modelexe"]
    else:
        modelexe = os.path.basename(o["model"]["data"]["modelexe"])
        model = o["model"]["data"]["modelexe"]
    try:
        os.makedirs(workdir)
    except:
        print(f"{workdir} exitsts")
    try:
        os.makedirs(outdir)
    except:
        print(f"{outdir} exitsts")
    try:
        os.makedirs(execdir)
    except:
        print(f"{execdir} exitsts")
    os.chdir(workdir)

    if not os.path.isdir(modeldir):
        print(f">>> Cannot find model directory {modeldir} <<<")
        sys.exit(1)

    if not os.path.isdir(utildir):
        print(f">>> Cannot find model directory {utildir} <<<")
        sys.exit(1)

    # Set srcdir based on modeldir
    srcdir = os.path.join(modeldir, 'src')

    # Check if srcdir exists
    if not os.path.isdir(srcdir):
        print(f">>> Cannot find model source directory {srcdir} <<<")
        sys.exit(1)

    # Convert srcdir to an absolute path
    srcdir = os.path.abspath(srcdir)  

    if tgcmdata == "None":
        tgcmdata = os.environ['TGCMDATA']
        print(f"Set tgcmdata = {tgcmdata}")

    if not os.path.isdir(tgcmdata):
        print(f">>> Cannot find data directory {tgcmdata}")

    # Check horizontal resolution
    if horires not in [5, 2.5, 1.25, 0.625]:
        print(f">>> Unknown model horizontal resolution {horires} <<<")
        sys.exit(1)
    
    # Check vertical resolution
    if vertres not in [0.5, 0.25, 0.125, 0.0625]:
        print(f">>> Unknown model vertical resolution {vertres} <<<")
        sys.exit(1)
    
    if nres_grid == "None" or nres_grid == None:
        if mres == 2:
            nres_grid = 5
        elif mres == 1:
            nres_grid = 6
        elif mres == 0.5:
            nres_grid = 7
        else:
            print(f">>> Unsupported magnetic resolution {mres} <<<")
            sys.exit(1)

    # Copy make if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, os.path.basename(make))):
        shutil.copy(os.path.join(utildir, os.path.basename(make)), execdir)
    # Copy Makefile if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, 'Makefile')):
        shutil.copy(os.path.join(utildir, 'Makefile'), execdir)

    # Copy mkdepends if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, 'mkdepends')):
        shutil.copy(os.path.join(utildir, 'mkdepends'), execdir)

    """
    if not os.path.isfile(input):
        print(f">>> Cannot find namelist input file {input} <<<")
        sys.exit(1)
    """
    
    if input == '' or output == '':
        input = os.path.abspath(input)
        output = os.path.abspath(output)
        
    util = os.path.abspath(utildir)


    coupling_file_path = os.path.join(execdir, 'coupling')

    # Check if the coupling file exists
    if os.path.isfile(coupling_file_path):
        with open(coupling_file_path, 'r') as file:
            lastcoupling = file.read().strip().lower() == 'true'

        # Compare coupling values
        if lastcoupling != coupling:
            print(f"Clean execdir {execdir} because coupling flag switched from {lastcoupling} to {coupling}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)
            with open(coupling_file_path, 'w') as file:
                file.write(str(coupling))
    else:
        # Create the coupling file and write the coupling value
        with open(coupling_file_path, 'w') as file:
            file.write(str(coupling))
        print(f"Created file coupling with coupling flag = {coupling}")


    debug_file_path = os.path.join(execdir, 'debug')

    # Check if the debug file exists
    if os.path.isfile(debug_file_path):
        with open(debug_file_path, 'r') as file:
            lastdebug = file.read().strip().lower() == 'true'

        # Compare debug values
        if lastdebug != debug:
            print(f"Clean execdir {execdir} because debug flag switched from {lastdebug} to {debug}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)

            with open(debug_file_path, 'w') as file:
                file.write(str(debug))
    else:
        # Create the debug file and write the debug value
        with open(debug_file_path, 'w') as file:
            file.write(str(debug))
        print(f"Created file debug with debug flag = {debug}")


    # Create the defs.h content
    defs_content = dedent(f"""\
    #define DLAT {horires}
    #define DLON {horires}
    #define GLON1 -180
    #define DLEV {vertres}
    #define ZIBOT -7
    #define ZITOP {zitop}
    #define NRES_GRID {nres_grid}
    """)

    # Write to defs.h
    defs_path = 'defs.h'
    with open(defs_path, 'w') as file:
        file.write(defs_content)

    # Check if defs.h exists in execdir and compare
    execdir_defs_path = os.path.join(execdir, 'defs.h')
    if os.path.isfile(execdir_defs_path):
        if not filecmp.cmp(defs_path, execdir_defs_path, shallow=False):
            # Files differ, switch resolutions
            print(f"Switching defs.h for model resolution {horires} x {vertres}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)
            shutil.copy(defs_path, execdir_defs_path)
        else:
            print(f"defs.h already set for model resolution {horires} x {vertres}")
    else:
        # defs.h does not exist in execdir, copy the file
        print(f"Copying {defs_path} to {execdir_defs_path} for resolution {horires} x {vertres}")
        shutil.copy(defs_path, execdir_defs_path)


    try:
        os.chdir(execdir)
        print(f"\nBegin building {model} in {os.getcwd()}")
    except OSError:
        print(f">>> Cannot cd to execdir {execdir}")
        sys.exit(1)



    # Create Make.env file
    make_env_path = os.path.join(execdir, 'Make.env')
    with open(make_env_path, 'w') as file:
        file.write(f"""MAKE_MACHINE  = {make}
    DIRS          = . {srcdir}
    EXECNAME      = {model}
    NAMELIST      = {input}
    OUTPUT        = {output}
    COUPLING      = {coupling}
    DEBUG         = {debug}
    """)

    # Build the model
    try:
        subprocess.run(['gmake', '-j8', 'all'], check=True)
        shutil.copy(model, workdir)
        print(f"Executable copied from {model} to {workdir}")
    except subprocess.CalledProcessError:
        print(">>> Error return from gmake all")
        sys.exit(1)

def create_pbs_scripts(options, run_name, segment_number):
    """
    Create PBS scripts for running TIEGCM model.

    Args:
        options (dict): A dictionary containing the model options.
        run_name (str): The name of the run.
        segment_number (int or None): The segment number of the run. If None, a single PBS script is created.

    Returns:
        str: The filepath of the created PBS script.

    Raises:
        FileNotFoundError: If the PBS template file is not found.

    """
    global PBS_TEMPLATE
    if PBS_TEMPLATE == None:
        PBS_TEMPLATE = os.path.join(options["model"]["data"]["modeldir"], 'tiegcmrun/template.pbs')
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)
    opt = copy.deepcopy(options)
    pbs_content = template.render(opt)
    workdir = opt["model"]["data"]["workdir"]
    if segment_number == None:
        pbs_script = os.path.join(workdir, f"{run_name}.pbs")
    else:
        pbs_script = os.path.join(workdir, f"{run_name}-{'{:02d}'.format(segment_number+1)}.pbs")
    with open(pbs_script, "w", encoding="utf-8") as f:
        f.write(pbs_content)
    return pbs_script

def create_inp_scripts(options, run_name, segment_number):
    """
    Create input scripts for running the TIEGCM model.

    Args:
        options (dict): A dictionary containing the model options.
        run_name (str): The name of the run.
        segment_number (int): The segment number.

    Returns:
        str: The path to the created input script.
    """
    global INP_TEMPLATE
    if INP_TEMPLATE == None:
        INP_TEMPLATE = os.path.join(options["model"]["data"]["modeldir"],'tiegcmrun/template.inp')
    with open(INP_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)
    opt = copy.deepcopy(options) 
    inp_content = template.render(opt)
    workdir = opt["model"]["data"]["workdir"]
    if segment_number == None:
        inp_script = os.path.join(workdir,f"{run_name}.inp")
    else:
        inp_script = os.path.join(workdir,f"{run_name}-{'{:02d}'.format(segment_number+1)}.inp")
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    with open(inp_script, "w", encoding="utf-8") as f:
        f.write(inp_content)    
    return inp_script

def find_file(pattern, path):
    """
    Find a file in the specified path that matches the given pattern. Assumes only one match.

    :param pattern: Pattern to look for in the file names.
    :param path: Path of the directory to search in.
    :return: File path if a match is found, else None.
    """
    for root, dirs, files in os.walk(path):  # Recursively go through all directories and subdirectories
        for name in files:
            if fnmatch.fnmatch(name, pattern):  # Check if file name matches the pattern
                return os.path.join(root, name)  # If so, return the file path immediately
    return None


def get_engage_start_time(datetime_str, seconds):
    # Parse the input datetime string
    dt = datetime.fromisoformat(datetime_str)
    
    # Subtract the seconds using timedelta
    new_dt = dt - timedelta(seconds=seconds)
    
    if new_dt.time() != datetime.min.time():
        # Adjust to the previous midnight
        new_dt = datetime.combine(new_dt.date(), datetime.min.time())
    
    # Return the new datetime as a string in the same format
    return new_dt.isoformat()

def time_to_dhms(time_str):

    # Convert string to datetime object
    time = datetime.strptime(time_str, "%Y-%m-%dT%H:%M:%S")
    
    # Extract day of year, hour of day, minute of hour, and second of minute
    day = time.timetuple().tm_yday
    hour = time.hour
    minute = time.minute
    second = time.second
    
    return [day, hour, minute, second]

def seconds_to_dhms(seconds):
    # Calculate the number of days
    days = seconds // (24 * 3600)
    seconds %= (24 * 3600)
    
    # Calculate the number of hours
    hours = seconds // 3600
    seconds %= 3600
    
    # Calculate the number of minutes
    minutes = seconds // 60
    
    # Calculate the remaining seconds
    seconds %= 60
    
    return [days, hours, minutes, seconds]

def gamres_to_res(gamres):
    "D", "Q", "O", "H"
    if gamres == "D":
        return 2.5 , 2.5
    elif gamres == "Q":
        return 2.5 , 1.25
    elif gamres == "O":
        return 2.5 , 0.625
    elif gamres == "H":
        return 2.5 ,0.625

def engage_parser(engage_parameters):
    """
    Parse the engage.json file and return the options dictionary.

    Args:
        engage_jsonfile (str): The path to the engage.json file.

    Returns:
        dict: The options dictionary.
    """

    o = engage_parameters["simulation"]
    
    hpc_system = o['hpc_system']
    coupled_job_name = o['job_name']
    coupled_start_date = o['start_date']
    stop_date = o['stop_date']

    use_segments = o['use_segments']
    segment_duration = int(float((o['segment_duration'])))
    segment = seconds_to_dhms(segment_duration)
    horires_standalone, horires_coupled = gamres_to_res(o['gamera_grid_type'])
    
    o = engage_parameters["pbs"]
    account_name = o['account_name']
    queue = o['queue']
    job_priority = o['job_priority']
    walltime = o['walltime']
    modules = o['modules']

    o = engage_parameters["coupling"]
    
    gamera_spin_up_time = int(o['gamera_spin_up_time'])
    gcm_spin_up_time = int(o['gcm_spin_up_time'])
    conda_env = o['conda_env']
    
    start_date = get_engage_start_time(coupled_start_date,gamera_spin_up_time+gcm_spin_up_time)
    
    o = engage_parameters["voltron"]
    voltron_dtOut = int(float(o['output']['dtOut']))
    hist = seconds_to_dhms(voltron_dtOut)
    STEP = int(float(o['coupling']['dtCouple']))

    root_directory= os.path.abspath(os.curdir)
    eo = engage_options = {}
    
    eo['job_name'] = coupled_job_name
    eo['hpc_system'] = hpc_system
    eo['start_time'] = start_date
    eo['coupled_start_time'] = coupled_start_date
    eo['stop_time'] = stop_date
    eo['segment'] = segment
    eo['segment_seconds'] = segment_duration
    eo['horires'] = horires_standalone
    eo['horires_coupled'] = horires_coupled
    eo['STEP'] = STEP
    eo["voltron_dtOut"] = voltron_dtOut
    eo['parentdir'] = root_directory

    eo['account_name'] = account_name
    eo['project_code'] = account_name
    eo['queue'] = queue
    eo['job_priority'] = job_priority
    eo['walltime'] = walltime
    eo['modules'] = modules
    eo['conda_env'] = conda_env

    eo['skip']= ['job_name','hpc_system','horires','parentdir','vertres', 'mres', 'input_file', 'LABEL','start_time','stop_time','secondary_start_time','secondary_stop_time','segment' ,'SOURCE_START','PRIHIST','MXHIST_PRIM','SECHIST','MXHIST_SECH','account_name','project_code','queue','job_priority','walltime']
    
    return engage_options


def segment_inp_pbs(options, run_name, pbs, engage_options=None):
    segment_times = segment_time(options["inp"]["start_time"], options["inp"]["stop_time"], [int(i) for i in options["inp"]["segment"].split()])
    pri_files = 0
    sec_files = 0
    og_options = copy.deepcopy(options)
    PRIHIST = og_options["inp"]["PRIHIST"]
    MXHIST_PRIM = og_options["inp"]["MXHIST_PRIM"]
    SECHIST = og_options["inp"]["SECHIST"]
    MXHIST_SECH = og_options["inp"]["MXHIST_SECH"]
    histdir = og_options["model"]["data"]["histdir"]
    workdir = og_options["model"]["data"]["workdir"]
    job_name = og_options["simulation"]["job_name"]
    inp_files = []
    pbs_files = []
    log_files = []
    pristart_times = []
    pristop_times = []
    with open(os.path.join(workdir,f'{run_name}.json'), "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)
    last_segment_time = len(segment_times) - 1 
    for segment_number, segment in enumerate(segment_times):
        segment_options = copy.deepcopy(og_options)
        segment_start = segment[0]
        segment_stop = segment[1]
        segment_options["simulation"]["job_name"] =  job_name+"-{:02d}".format(segment_number +1)
        if segment_number == 0:
            segment_options["inp"]["SOURCE"] = og_options["inp"]["SOURCE"]
            segment_options["inp"]["SOURCE_START"] = og_options["inp"]["SOURCE_START"]
            segment_START_YEAR, segment_START_DAY, segment_PRISTART, segment_PRISTOP = inp_pri_date(segment_start,segment_stop)
            segment_options["inp"]["START_YEAR"] = segment_START_YEAR
            segment_options["inp"]["START_DAY"] = segment_START_DAY
            segment_options["inp"]["PRISTART"] = ' '.join(map(str, segment_PRISTART))
            segment_options["inp"]["PRISTOP"] = ' '.join(map(str, segment_PRISTOP))
            segment_OUTPUT, pri_files = inp_pri_out(segment_start, segment_stop, [int(i) for i in PRIHIST.split()], MXHIST_PRIM, pri_files, histdir,run_name)
            segment_options["inp"]["OUTPUT"] = segment_OUTPUT
        else:
            segment_START_YEAR, segment_START_DAY, segment_PRISTART, segment_PRISTOP = inp_pri_date(segment_start,segment_stop)
            segment_options["inp"]["SOURCE"] = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
            segment_options["inp"]["SOURCE_START"] = ' '.join(map(str, segment_PRISTART))
            segment_options["inp"]["START_YEAR"] = segment_START_YEAR
            segment_options["inp"]["START_DAY"] = segment_START_DAY
            segment_options["inp"]["PRISTART"] = ' '.join(map(str, segment_PRISTART))
            segment_options["inp"]["PRISTOP"] = ' '.join(map(str, segment_PRISTOP))
            if segment_number == last_segment_time and engage_options != None:
                segment_START_YEAR, segment_START_DAY, segment_PRISTART, segment_PRISTOP = inp_pri_date(segment_start,segment_stop)
                segment_PRISTOP_day, segment_PRISTOP_hour, segment_PRISTOP_minute, segment_PRISTOP_second = segment_PRISTOP
                if segment_PRISTOP_second != 0:
                    segment_PRIHIST = [0, 0, 0, 1]
                    segment_MXHIST_PRIM = segment_PRISTOP_minute*60 + segment_PRISTOP_hour*60 + segment_PRISTOP_second
                elif segment_PRISTOP_minute != 0:
                    segment_PRIHIST = [0, 0, 1, 0]
                    segment_MXHIST_PRIM = segment_PRISTOP_hour*60 + segment_PRISTOP_minute
                elif segment_PRISTOP_hour != 0:
                    segment_PRIHIST = [0, 1, 0, 0]
                    segment_MXHIST_PRIM = segment_PRISTOP_hour
                else:
                    segment_PRIHIST = PRIHIST
                    segment_MXHIST_PRIM = MXHIST_PRIM
                segment_OUTPUT, pri_files = inp_pri_out(segment_start, segment_stop, segment_PRIHIST, segment_MXHIST_PRIM, pri_files, histdir,run_name)
                segment_options["inp"]["PRIHIST"] = ' '.join(map(str, segment_PRIHIST))
                segment_options["inp"]["MXHIST_PRIM"] = segment_MXHIST_PRIM
                segment_options["inp"]["OUTPUT"] = segment_OUTPUT
            else:
                segment_OUTPUT, pri_files = inp_pri_out(segment_start, segment_stop, [int(i) for i in PRIHIST.split()], MXHIST_PRIM, pri_files, histdir,run_name)
                segment_options["inp"]["OUTPUT"] = segment_OUTPUT
        segment_SECSTART, segment_SECSTOP = inp_sec_date(segment_start, segment_stop, [int(i) for i in SECHIST.split()])
        segment_options["inp"]["SECSTART"] = ' '.join(map(str, segment_SECSTART))
        segment_options["inp"]["SECSTOP"] = ' '.join(map(str, segment_SECSTOP))
        segment_SECOUT, sec_files = inp_sec_out(segment_start, segment_stop, [int(i) for i in SECHIST.split()], MXHIST_SECH, sec_files, histdir,run_name)
        segment_options["inp"]["SECOUT"] = segment_SECOUT
        segment_options["model"]["data"]["input_file"] = create_inp_scripts(segment_options,run_name,segment_number)
        pristart_times.append(segment_options["inp"]["PRISTART"])
        pristop_times.append(segment_options["inp"]["PRISTOP"])
        if segment_number == 0:
            init_inp = segment_options["model"]["data"]["input_file"]
        segment_options["model"]["data"]["log_file"] = os.path.join( options["model"]["data"]["workdir"],f"{run_name}-{'{:02d}'.format(segment_number+1)}.out")
        if pbs == True:
            if options["simulation"]["hpc_system"] != "linux":
                '''
                if segment_number != len(segment_times) - 1:
                    next_pbs = os.path.join(workdir,f"{run_name}-{'{:02d}'.format(segment_number+1)}.pbs")
                    segment_options["job"]["job_chain"] = [f"qsub {next_pbs}"]
                    segment_options["job"]["job_chain"] = [None]
                else:
                    segment_options["job"]["job_chain"] = [None]
                '''
                if segment_number == last_segment_time and engage_options != None:
                    interpolation_script = os.path.join(segment_options["model"]["data"]["workdir"],f'interpolation.py')
                    with open(interpolation_script, "w", encoding="utf-8") as f:
                        f.write("import sys\n")
                        f.write(f"sys.path.append('{TIEGCMHOME}/tiegcmrun')\n")
                        f.write("import tiegcmrun\n")
                        f.write("print(f'tiegcmrum from {tiegcmrun.__file__}')\n")
                        horires_coupled = engage_options["horires_coupled"]
                        vertres_coupled, mres_coupled, nres_grid_coupled, STEP_coupled = resolution_solver(horires_coupled,engage_options)
                        SOURCE_coupling = os.path.join(os.path.dirname(segment_options["model"]["data"]["workdir"]),f'{engage_options["job_name"]}_prim.nc')
                        input_standalone = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
                        f.write(f"tiegcmrun.interpic('{input_standalone}',{float(horires_coupled)},{float(vertres_coupled)},{float(segment_options['model']['specification']['zitop'])},'{SOURCE_coupling}')\n")
                        interpolation_pbs = [f'conda activate {engage_options["conda_env"]}',f'python {interpolation_script}']
                    segment_options["job"]["job_chain"] = interpolation_pbs
                pbs_script = create_pbs_scripts(segment_options,run_name,segment_number)
                if segment_number == 0:
                    init_pbs = pbs_script
            else:
                pbs_script = None
        else:
            pbs_script = None
        inp_files.append(segment_options["model"]["data"]["input_file"])
        pbs_files.append(pbs_script)
        log_files.append(segment_options["model"]["data"]["log_file"])
    return inp_files, pbs_files,log_files, pristart_times, pristop_times

def engage_run(options, debug, coupling, engage):
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)
    options_standalone = copy.deepcopy(options)
    options_coupling = copy.deepcopy(options)
    #For standalone
    pbs=True
    options_standalone["simulation"]["job_name"] = f'{engage["job_name"]}-tiegcm-standalone'
    options_standalone["inp"]["stop_time"] = engage["coupled_start_time"]
    options_standalone["inp"]["PRIHIST"] = '1 0 0 0'
    options_standalone["inp"]["MXHIST_PRIM"] = 1
    options_standalone["inp"]["SECHIST"] = '0 1 0 0'
    options_standalone["inp"]["MXHIST_SECH"] = 24
    options_standalone["inp"]["segment"] = '1 0 0 0'
    options_standalone["model"]["data"]["workdir"] = os.path.join(engage["parentdir"],"tiegcm_standalone")
    if not os.path.exists(options_standalone["model"]["data"]["workdir"]):
        os.makedirs(options_standalone["model"]["data"]["workdir"])
    options_standalone["model"]["data"]["histdir"] = os.path.join(engage["parentdir"],"tiegcm_standalone")
    
    in_prim = options_standalone["inp"]["SOURCE"]
    out_prim = f'{options_standalone["model"]["data"]["workdir"]}/{options_standalone["simulation"]["job_name"]}_prim.nc'
    options_standalone["inp"]["SOURCE"] = out_prim
    horires_standalone= engage["horires"]
    vertres_standalone, mres_standalone, nres_grid_standalone, STEP_standalone = resolution_solver(horires_standalone,engage)
    interpic (in_prim,float(horires_standalone),float(vertres_standalone),float(options_standalone['model']['specification']['zitop']),out_prim)    
    standalone_inp_files,standalone_pbs_files, standalone_log_files,pristart_times, pristop_times=segment_inp_pbs(options_standalone, options_standalone["simulation"]["job_name"],pbs, engage)
    #For coupled
    pbs=False
    options_coupling["model"]["data"]["modelexe"] = options_coupling["model"]["data"]["coupled_modelexe"]
    coupling_modelexe=options_coupling["model"]["data"]["coupled_modelexe"]
    horires_coupling = float(engage["horires_coupled"])
    options_coupling["model"]["specification"]["horires"] = horires_coupling
    vertres_coupling, mres_coupling, nres_grid_coupling, STEP_coupling = resolution_solver(horires_coupling,engage)
    options_coupling["model"]["specification"]["vertres"] = vertres_coupling
    options_coupling["model"]["specification"]["mres"] = mres_coupling
    options_coupling["model"]["specification"]["nres_grid"] = nres_grid_coupling
    options_coupling["inp"]["STEP"] = STEP_coupling
    SOURCE_coupling = os.path.join(options_coupling["model"]["data"]["workdir"],f'{engage["job_name"]}_prim.nc')
    options_coupling["inp"]["SOURCE"] = SOURCE_coupling #standalone_inp_files[-1]
    options_coupling["inp"]["SOURCE_START"] = pristop_times[-1]
    options_coupling["simulation"]["job_name"] = f'{engage["job_name"]}'
    options_coupling["inp"]["start_time"] = engage["coupled_start_time"]
    options_coupling["inp"]["PRIHIST"] = " ".join(str(i) for i in engage["segment"])
    options_coupling["inp"]["MXHIST_PRIM"] = 1
    options_coupling["inp"]["SECHIST"] = " ".join(str(i) for i in seconds_to_dhms(engage["voltron_dtOut"]))
    options_coupling["inp"]["MXHIST_SECH"] = int(engage["segment_seconds"]/engage["voltron_dtOut"])
    options_coupling["inp"]["GSWM_MI_DI_NCFILE"] = find_file(f'*gswm_diurn_{horires_coupling}d_99km*', TIEGCMDATA)
    options_coupling["inp"]["GSWM_MI_SDI_NCFILE"] = find_file(f'*gswm_semi_{horires_coupling}d_99km*', TIEGCMDATA)
    options_coupling["inp"]["GSWM_NM_DI_NCFILE"] = find_file(f'*gswm_nonmig_diurn_{horires_coupling}d_99km*', TIEGCMDATA)
    options_coupling["inp"]["GSWM_NM_SDI_NCFILE"] = find_file(f'*gswm_nonmig_semi_{horires_coupling}d_99km*', TIEGCMDATA)
    coupling_inp_files,coupling_pbs_files, coupling_log_files, pristart_times, pristop_times = segment_inp_pbs(options_coupling, options_coupling["simulation"]["job_name"],pbs)
    select_coupling,ncpus_coupling,mpiprocs_coupling=select_resource_defaults(options_coupling,option_descriptions)
    options_coupling["job"]["resource"]["select"] = select_coupling
    options_coupling["job"]["resource"]["ncpus"] = ncpus_coupling
    options_coupling["job"]["resource"]["mpiprocs"] = mpiprocs_coupling
    nprocs_coupling = int(mpiprocs_coupling)*int(select_coupling)
    options_coupling["job"]["nprocs"] = nprocs_coupling
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)
    
    return options_coupling,standalone_pbs_files,coupling_inp_files

def tiegcmrun(args=None):
    # Set up the command-line parser.
    parser = create_command_line_parser()
    if args is not None:
        args = parser.parse_args(args)
    else:
        args = parser.parse_args()
    clobber = args.clobber
    debug = args.debug
    options_path = args.options_path
    verbose = args.verbose
    coupling = args.coupling
    compile = args.compile
    onlycompile = args.onlycompile
    execute = args.execute
    benchmark = args.benchmark
    engage = args.engage
    if args.engage != None:
        args.engage = engage_parser(json.loads(engage))
    mode = args.mode   
    # Fetch the run options.
    if benchmark != None and mode == None:
        mode = "BENCH"
    elif mode == None:
        mode = "BASIC"
    if compile == True or onlycompile == True:
        compile_flag = True
    elif compile == False and onlycompile == False:
        compile_flag = False
    print("\n")
    print("Instructions:")
    print(f"-> Default Selected input parameter is given in {GREEN}GREEN{RESET}")
    print(f"-> Warnings and Information are given in {YELLOW}YELLOW{RESET}")
    print(f"-> Errors are given in {RED}RED{RESET}")
    print(f"-> Valid values (if any) are given in brackets eg. (value1 | value2 | value3) ")
    print(f"-> Enter '?' for any input parameter to get a detailed description")
    print(f"\n")
    print("Run Options:")
    if benchmark != None:
        print(f"Benchmark = {benchmark}")
    print(f"User Mode = {mode}")
    print(f"Compile = {compile_flag}")
    print(f"Execute = {execute}")
    print(f"Coupling = {coupling}")  
    if args.engage != None:
        print(f"Engage = True")
    print(f"\n") 

    if options_path:
        # Read the run options from a JSON file.
        with open(options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
    else:
        # Prompt the user for the run options.
        options = prompt_user_for_run_options(args)
    if debug:
        print(f"options = {options}")

    # Move to the run directory.
    run_name = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
    execdir   = options["model"]["data"]["execdir"]
    workdir = options["model"]["data"]["workdir"]
    outdir = options["model"]["data"]["histdir"]
    try:
        os.makedirs(workdir)
    except:
        print(f"{workdir} exitsts")
    try:
        os.makedirs(outdir)
    except:
        print(f"{outdir} exitsts")
    try:
        os.makedirs(execdir)
    except:
        print(f"{execdir} exitsts")
    # Save the options dictionary as a JSON file in the current directory.
    path = f"{workdir}/{run_name}.json"
    
    tiegcmdata = options["model"]["data"]["tgcmdata"]
    horires = options["model"]["specification"]["horires"]
    vertres = options["model"]["specification"]["vertres"]
    zitop = options["model"]["specification"]["zitop"]
    if options.get("inp") == None:
        input_file_generatred = True
    else:
        input_file_generatred = False
    if args.onlycompile == True:
        compile_tiegcm(options, debug, coupling)
    elif args.engage != None:
        """
        if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
            in_prim = options["inp"]["SOURCE"]
            out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
            options["inp"]["SOURCE"] = out_prim
            interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
        else:
            out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
            options["inp"]["SOURCE"] = out_prim
            print(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc exists')
        """
        options_coupling,standalone_pbs_files,coupling_inp_files = engage_run(options, debug, coupling, args.engage)
        return (options_coupling,standalone_pbs_files,coupling_inp_files)
    else:
        if options["model"]["specification"]["segmentation"] == False:
            if options["model"]["data"]["input_file"] == None or not os.path.isfile(options["model"]["data"]["input_file"]):
                if input_file_generatred == False:
                    if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
                        if args.benchmark == None:
                            in_prim = options["inp"]["SOURCE"]
                            out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                            options["inp"]["SOURCE"] = out_prim
                            interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
                        else:
                            in_prim = options["inp"]["SOURCE"]
                            out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                            options["inp"]["SOURCE"] = out_prim
                            interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
                    else:
                        out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                        options["inp"]["SOURCE"] = out_prim
                        print(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc exists')
                    options["model"]["data"]["input_file"] = create_inp_scripts(options,run_name,None)
            if options["model"]["data"]["log_file"] == None:
                options["model"]["data"]["log_file"] = os.path.join( options["model"]["data"]["workdir"], f"{run_name}.out")

            if os.path.exists(path):
                if not clobber:
                    raise FileExistsError(f"Options file {path} exists!")
            with open(path, "w", encoding="utf-8") as f:
                json.dump(options, f, indent=JSON_INDENT)

            if args.compile == True:
                compile_tiegcm(options, debug, coupling)
            
            if options["simulation"]["hpc_system"] != "linux":
                pbs_script = create_pbs_scripts(options,run_name,None)
            if args.execute == True and options["simulation"]["hpc_system"] != "linux":
                if args.compile == False:
                    if find_file(options["model"]["data"]["modelexe"], execdir) == None and os.path.exists(options["model"]["data"]["modelexe"]) == False :
                        print(f'{RED}Unable to find executable in {execdir}{RESET}')
                        exit(1)
                try:
                    result = subprocess.run(['qsub', pbs_script], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    job_id = result.stdout.strip()
                    print(f'Job submitted successfully. Job ID: {job_id}')
                except subprocess.CalledProcessError as e:
                    print(f'{YELLOW}Error submitting job: {e.stderr}{RESET}')
                    print(f"{YELLOW}Check PBS script for erros{RESET}")
                    print(f"To submit job use command {YELLOW}qsub {pbs_script}{RESET}")
            elif options["simulation"]["hpc_system"] != "linux":
                print(f"{YELLOW}Execute is set to false{RESET}")
                print(f"{YELLOW}To submit job use command{RESET} qsub {pbs_script}")
            else:
                print(f"{YELLOW}HPC System is set to linux{RESET}")
                print(f"{YELLOW}To run the model use command{RESET} mpirun {options['model']['data']['modelexe']} {options['model']['data']['input_file']}")
        else:
            if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
                if args.benchmark == None:
                    in_prim = options["inp"]["SOURCE"]
                    out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                    options["inp"]["SOURCE"] = out_prim
                    interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
                else:
                    in_prim = options["inp"]["SOURCE"]
                    out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                    options["inp"]["SOURCE"] = out_prim
                    interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
            else:
                out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                options["inp"]["SOURCE"] = out_prim
                print(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc exists')
            inp_files, pbs_files, pristart_times = segment_inp_pbs(options, run_name, pbs=True)
            init_inp = inp_files[0]    
            init_pbs = pbs_files[0]
            options["model"]["data"]["input_file"] = init_inp
            if args.compile == True:
                compile_tiegcm(options, debug, coupling)

            if args.execute == True:
                if args.compile == False:
                    if find_file(options["model"]["data"]["modelexe"], execdir) == None:
                        print(f'{RED}Unable to find executable in {execdir}{RESET}')
                        exit(1)
                try:
                    result = subprocess.run(['qsub', init_pbs], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    job_id = result.stdout.strip()
                    print(f'Job submitted successfully. Job ID: {job_id}')
                except subprocess.CalledProcessError as e:
                    print(f'{YELLOW}Error submitting job: {e.stderr}{RESET}')
                    print(f"{YELLOW}Check PBS script for erros{RESET}")
                    print(f"To submit job use command {YELLOW}qsub {pbs_script}{RESET}")
            else:
                print(f"{YELLOW}Execute is set to false{RESET}")
                print(f"{YELLOW}To submit job use command{RESET} qsub {pbs_script}")

if __name__ == "__main__":
    """Begin tiegcmrun program."""
    tiegcmrun()
