#!/usr/bin/env python


"""makeitso for the MAGE magnetosphere software.

This script is perforsm all of the steps needed to prepare to run a MAGE
magnetosphere simulation run. By default, this script is interactive - the user
is prompted for each decision  that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BASIC" (the default) - the user is prompted to set only a small subset of MAGE
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.
"""


# Import standard modules.
import argparse
import copy
import datetime
import json
import os
import subprocess
import sys
import shutil
import subprocess
import filecmp
from textwrap import dedent
from numpy import ndarray, interp, log, exp, linspace, allclose, mean
from netCDF4 import Dataset
from argparse import ArgumentParser
from os.path import isfile, splitext
import fnmatch
# Import 3rd-party modules.

from jinja2 import Template

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare a MAGE magnetosphere model run."

# Indent level for JSON output.
JSON_INDENT = 4

# Path to current tiegcm datafiles
try:
    TIEGCMDATA = os.environ["TIEGCMDATA"]
except:
    os.environ['TIEGCMDATA'] = input(f'Enter TIEGCM data directory [TIEGCMDATA]: ')
    TIEGCMDATA = os.environ["TIEGCMDATA"]

# Path to current tiegcm installation
try:
    TIEGCMHOME = os.environ["TIEGCMHOME"]
except:
    os.environ['TIEGCMHOME'] = input(f'Enter TIEGCM model directory [TIEGCMHOME]: ')
    TIEGCMHOME = os.environ["TIEGCMHOME"]
# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(TIEGCMHOME, "tiegcmrun")

# Path to template .inp file.
INP_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.inp")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")

BENCHAMRKS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, 'benchmarks.json')


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
    lower_cap = 1e-8
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
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default="BASIC",
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
        "--execute","-e", action="store_false",
        help="Execute TIEGCM (default: %(default)s)."
    )
    parser.add_argument(
        "--compile","-c", action="store_false",
        help="Compile TIEGCM (default: %(default)s)."
    )
    parser.add_argument(
        "--benchmark", "-bench", default = None, type=valid_bench,
        help="Bechmark run name  (default: %(default)s)."
    )
    return parser

def valid_bench(value):
    # Custom validation logic
    if value not in [ None, 
                     'seasons', 'decsol_smax', 'decsol_smin', 'junsol_smax', 'junsol_smin','mareqx_smax', 'mareqx_smin', 'sepeqx_smax', 'sepeqx_smin',
                     'storms', 'dec2006_heelis_gpi', 'dec2006_weimer_imf', 'jul2000_heelis_gpi', 'jul2000_weimer_imf', 'nov2003_heelis_gpi', 'nov2003_weimer_imf', 'whi2008_heelis_gpi', 'whi2008_weimer_imf',
                     'climatology', 'climatology_smax', 'climatology_smin' 
                    ]:
        raise argparse.ArgumentTypeError(f"{value} is not a valid benchmark option.")
    return value

def get_run_option(name, description, mode="BASIC"):
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
    # Compare the current mode to the parameter level setting. If the variable
    # level is higher than the user mode, just use the default.
    if mode == "NOVICE" and level in ["BASIC","INTERMEDIATE", "EXPERT"]:
        return default
    if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
        return default
    if mode == "INTERMEDIATE" and level in ["EXPERT"]:
        return default

    # If provided, add the valid values in val1|val2 format to the prompt.
    if valids is not None:
        vs = "|".join(valids)
        prompt += f" ({vs})"

    # If provided, add the default to the prompt.
    if default is not None:
        prompt += f" [{default}]"

    # Prompt the user and fetch the input until a good value is provided.
    ok = False
    option_value = ""
    first_pass = True
    while not ok:
        if name in ("SECFLDS","GSWM_data" ,"other_input", "other_pbs","other_job"):
            temp_value = input(f"{prompt} or enter: ")
            if temp_value != "":
                option_value = option_value +'"' +temp_value + '",'
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
            else:
                if option_value != "":
                    option_value = "["+ option_value[:-1] + "]"
                    option_value = json.loads(option_value)
                    ok = True
        else:
            # Fetch input from the user.
            option_value = input(f"{prompt}: ")

            # Use the default if no user input provided.
            if option_value == "":
                option_value = default
            # Use None if user input is none or None.
            elif option_value == 'none' or option_value == 'None':
                option_value = None
            # Validate the result. If bad, start over.
            if valids is not None and option_value not in valids:
                print(f"Invalid value for option {name}: {option_value}!")
                continue

            # Keep this result.
            ok = True
            if option_value != None:
                option_value = str(option_value)
    # Return the option as a string.
    return option_value


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
    # Save the user mode.
    mode = args.mode
    benchmark = args.benchmark
    if benchmark != None and mode == 'BASIC':
        mode = "NOVICE"
    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)
    with open(BENCHAMRKS_FILE, "r", encoding="utf-8") as f:
        benchmarks_options = json.load(f)
    # Initialize the dictionary of program options.
    options = {}

    #-------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_descriptions["simulation"]
    if benchmark != None:
        od["job_name"]["default"] = benchmark
    # Prompt for the parameters.
    for on in ["job_name", "hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)
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
    o["modeldir"] = get_run_option("modeldir", od["modeldir"], mode)
    od["parentdir"]["default"] = os.getcwd() 
    o["parentdir"] = get_run_option("parentdir", od["parentdir"], mode)
    if o["parentdir"] == None:
        temp_mode = "INTERMEDIATE"
        o["execdir"] = get_run_option("execdir", od["execdir"], temp_mode)
        od["workdir"]["default"] = o["execdir"]
        od["histdir"]["default"] = o["execdir"]
    else:
        od["execdir"]["default"] = os.path.join(o["parentdir"],"exec")
        od["workdir"]["default"] = os.path.join(o["parentdir"],"stdout")
        od["histdir"]["default"] = os.path.join(o["parentdir"],"hist")
        o["execdir"] = get_run_option("execdir", od["execdir"], mode)

    o["workdir"] = get_run_option("workdir", od["workdir"], temp_mode)
    o["histdir"] = get_run_option("histdir", od["histdir"], temp_mode)
    temp_mode = mode    
    o["utildir"] = os.path.join(o["modeldir"],'scripts')
    od["tgcmdata"]["default"] = TIEGCMDATA
    o["tgcmdata"] = get_run_option("tgcmdata", od["tgcmdata"], mode)
    o["input_file"] = get_run_option("input_file", od["input_file"], mode)
    o["output_file"] = get_run_option("output_file", od["output_file"], mode)

    if od["make"]["default"] == None:
        if options["simulation"]["hpc_system"] == "derecho":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_de')
        elif options["simulation"]["hpc_system"] == "pleiades":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_pf')
        elif options["simulation"]["hpc_system"] == "None":
            od["make"]["default"] = os.path.join(options["model"]["data"]["modeldir"],'scripts/Make.intel_linux')
    o["make"] = get_run_option("make", od["make"], mode)

    o["modelexe"] = get_run_option("modelexe", od["modelexe"], mode)

    tiegcmdata_dir = o["tgcmdata"]
    #------------------------------------

    # Specification Options
            
    options["model"]["specification"] = {}
    o = options["model"]["specification"]

    od = option_descriptions["model"]["specification"]                
    for on in od:
        if on == "vertres":
            if float(horires) == 5:
                od["vertres"]["default"] = 0.5
            elif float(horires) == 2.5:
                od["vertres"]["default"] = 0.25
            elif float(horires) == 1.25:
                od["vertres"]["default"] = 0.125
            elif float(horires) == 0.625:
                od["vertres"]["default"] = 0.0625
        o[on] = get_run_option(on, od[on], mode)
        if on =="horires":
            horires = float(o[on])

    if o["nres_grid"] == None:
        if float(o["mres"]) == 2:
            o["nres_grid"] = 5
        elif float(o["mres"]) == 1:
            o["nres_grid"] = 6
        elif float(o["mres"]) == 0.5:
            o["nres_grid"] = 7

    #-------------------------------------------------------------------------
    # INP options
    options["inp"] = {}
    o = options["inp"]

    od = option_descriptions["inp"]     
    

    if float(horires) == 5:
        od["STEP"]["default"] = 60
    elif float(horires) == 2.5:
        od["STEP"]["default"] = 30
    elif float(horires) == 1.25:
        od["STEP"]["default"] = 10
    elif float(horires) == 0.625:
        od["STEP"]["default"] = 5
    run_name = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
    histdir = options["model"]["data"]["histdir"]
    if benchmark != None:
        oben = benchmarks_options[benchmark]["inp"]
        for on in od:
            if on in oben:
                if on in ["OUTPUT", "SECOUT"]:
                    temp_output = oben[on]
                    temp_output = temp_output.replace("+histdir+", histdir)
                    temp_output = temp_output.replace("+run_name+", run_name)
                    od[on]["default"] = temp_output
                elif on in ["GSWM_data", "other_input"]:
                        temp_output = [item.replace("+tiegcmdata+", tiegcmdata_dir) if item is not None and item != 'null' else item for item in oben[on]]
                        od[on]["default"] = oben[on]
                elif oben[on] == None:
                    od[on]["default"]
                else:
                    od[on]["default"] = oben[on]

    od["LABEL"]["default"] = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"

    temp_mode = mode

    for on in od:
        if on == "start_date":
            o[on] = get_run_option(on, od[on], temp_mode)
            if o[on] == None and benchmark == None:
                temp_mode =  "INTERMEDIATE"
        elif on == "SECSTART":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "SECSTART = "+str(temp_val)
        elif on == "SECSTOP":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "SECSTOP = "+str(temp_val)
        elif on == "SECHIST":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "SECHIST = "+str(temp_val)
        elif on == "SECOUT":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "SECOUT = "+str(temp_val)
        elif on == "MXHIST_SECH":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "MXHIST_SECH = "+str(temp_val)
        elif on == "SECFLDS":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = temp_val
        elif on == "POWER":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "POWER = "+str(temp_val)
        elif on == "CTPOTEN":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "CTPOTEN = "+str(temp_val)
        elif on == "F107":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "F107 = "+str(temp_val)
        elif on == "F107A":
            temp_val = get_run_option(on, od[on], temp_mode)
            if temp_val != None:
                o[on] = "F107A = "+str(temp_val)    
        elif on == "GSWM_data":
            o[on] = get_run_option(on, od[on], temp_mode)
            if o[on] == [None]:
                if  float(options["model"]["specification"]["horires"]) == 2.5:
                    o[on] = [
                        f"GSWM_MI_DI_NCFILE = '{tiegcmdata_dir}/gswm_diurn_2.5d_99km.nc'",
                        f"GSWM_MI_SDI_NCFILE = '{tiegcmdata_dir}/gswm_semi_2.5d_99km.nc'",
                        f"GSWM_NM_DI_NCFILE = '{tiegcmdata_dir}/gswm_nonmig_diurn_2.5d_99km.nc'",
                        f"GSWM_NM_SDI_NCFILE = '{tiegcmdata_dir}/gswm_nonmig_semi_2.5d_99km.nc'"
                    ]
        else:
            o[on] = get_run_option(on, od[on], temp_mode)
            

    #-------------------------------------------------------------------------

    # PBS options
    options["pbs"] = {}
    o = options["pbs"]

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # HPC platform-specific options
    hpc_platform = options["simulation"]["hpc_system"]
    od = option_descriptions["pbs"][hpc_platform]
    if hpc_platform == "derecho":
        odt=od["resource"]
        for on in odt:
            if on == "select":
                if float(horires) == 5:
                    odt["select"]["default"] = 1
                elif float(horires) == 2.5:
                    odt["select"]["default"] = 1
                elif float(horires) == 1.25:
                    odt["select"]["default"] = 2
                elif float(horires) == 0.625:
                    odt["select"]["default"] = 3
            if on == "ncpus":
                if float(horires) == 5:
                    odt["ncpus"]["default"] = 72
                elif float(horires) == 2.5:
                    odt["ncpus"]["default"] = 72
                elif float(horires) == 1.25:
                    odt["ncpus"]["default"] = 72
                elif float(horires) == 0.625:
                    odt["ncpus"]["default"] = 96
    for on in od:
        if on == "resource":
            options["pbs"]["resource"] = {}
            odt = od["resource"]
            ot = options["pbs"]["resource"]
            for ont in odt:

                if ont == "mpiprocs":
                    ot[ont] = get_run_option(ont, odt[ont], mode)
                    if ot[ont] == "ncpus":
                        ot[ont] = ot["ncpus"]
                    mpiprocs = ot[ont]
                ot[ont] = get_run_option(ont, odt[ont], mode)
                if ont == "select":
                    nnodes = ot[ont]
                elif ont == "ncpus":
                    ncpus = ot[ont]
                ot[ont] = ont+"="+str(ot[ont])
        elif on == "job_priority":
            o[on] = "#PBS -l job_priority="+get_run_option(on, od[on], mode)
        else:
            o[on] = get_run_option(on, od[on], mode)

    o["nprocs"] = int(nnodes) * int(mpiprocs)
    # Return the options dictionary.
    return options

def compile_tiegcm(options, debug, coupling):
    o = options
    modeldir  = o["model"]["data"]["modeldir"]
    execdir   = o["model"]["data"]["execdir"]
    workdir = o["model"]["data"]["workdir"]
    outdir = o["model"]["data"]["histdir"]
    tgcmdata  = o["model"]["data"]["tgcmdata"]
    utildir   = os.path.join(o["model"]["data"]["modeldir"],"scripts")
    input     = o["model"]["data"]["input_file"]
    output    = o["model"]["data"]["output_file"]
    horires   = float(o["model"]["specification"]["horires"])
    vertres   = float(o["model"]["specification"]["vertres"])
    zitop     = float(o["model"]["specification"]["zitop"])
    mres      = float(o["model"]["specification"]["mres"])
    nres_grid = float(o["model"]["specification"]["nres_grid"])
    make      = o["model"]["data"]["make"]
    coupling  = False
    modelexe = o["model"]["data"]["modelexe"]

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

    if not os.path.isfile(input):
        print(f">>> Cannot find namelist input file {input} <<<")
        sys.exit(1)
    
    model = os.path.join(execdir, modelexe)
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
    except subprocess.CalledProcessError:
        print(">>> Error return from gmake all")
        sys.exit(1)

def create_pbs_scripts(options):
    global PBS_TEMPLATE
    if PBS_TEMPLATE == None:
        PBS_TEMPLATE = os.path.join(options["model"]["data"]["modeldir"],'tiegcmrun/template.pbs')
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)
    opt = copy.deepcopy(options) 
    pbs_content = template.render(opt)
    pbs_script = os.path.join(
            opt["model"]["data"]["workdir"],
            f"{opt['simulation']['job_name']}_{opt['model']['specification']['horires']}x{opt['model']['specification']['vertres']}.pbs"
        )
    with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)
    return pbs_script

def create_inp_scripts(options):
    global INP_TEMPLATE
    if INP_TEMPLATE == None:
        INP_TEMPLATE = os.path.join(options["model"]["data"]["modeldir"],'tiegcmrun/template.inp')
    with open(INP_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)
    opt = copy.deepcopy(options) 
    inp_content = template.render(opt)
    workdir = opt["model"]["data"]["workdir"]
    inp_script = os.path.join(
            workdir,
            f"{opt['simulation']['job_name']}_{opt['model']['specification']['horires']}x{opt['model']['specification']['vertres']}.inp"
        )
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
def main():
    # Set up the command-line parser.
    parser = create_command_line_parser()
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    clobber = args.clobber
    debug = args.debug
    options_path = args.options_path
    verbose = args.verbose
    coupling = args.coupling    
    # Fetch the run options.
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
    #os.chdir(options["pbs"]["run_directory"])
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
    if options["model"]["data"]["input_file"] == None or not os.path.isfile(options["model"]["data"]["input_file"]):
        if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
            if args.benchmark == None:
                in_prim = options["inp"]["SOURCE"]
                out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                options["inp"]["SOURCE"] = out_prim
                interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
            else:
                prim_loc = os.path.join(tiegcmdata,'prim')
                in_prim = find_file('*'+options['simulation']['job_name']+'*', prim_loc)
                out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
                options["inp"]["SOURCE"] = out_prim
        else:
            print(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc exists')
        options["model"]["data"]["input_file"] = create_inp_scripts(options)
    if options["model"]["data"]["output_file"] == None:
        options["model"]["data"]["output_file"] = os.path.join( options["model"]["data"]["workdir"], f"{run_name}.out")


    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

    if args.compile == True:
        compile_tiegcm(options, debug, coupling)
    
    pbs_script = create_pbs_scripts(options)

    if args.execute == True:
        try:
            result = subprocess.run(['qsub', pbs_script], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            job_id = result.stdout.strip()
            print(f'Job submitted successfully. Job ID: {job_id}')
        except subprocess.CalledProcessError as e:
            print(f'Error submitting job: {e.stderr}')


if __name__ == "__main__":
    """Begin main program."""
    main()