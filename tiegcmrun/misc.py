"""
misc.py for the TIEGCMrun software.

This script contains various utility functions that are used to support the TIEGCMrun software. These functions include file handling, time segmentation, validation, resolution solving, and default selection for various options.

Functions included:

- get_mtime(file_path): Extracts and pads 'mtime' data from a given file.
- segment_time(start_time_str, stop_time_str, interval_array): Generates a list of time intervals between a start and stop time.
- valid_bench(value): Validates the benchmark option.
- resolution_solver(horires, engage_options=None): Determines vertical resolution, model resolution, grid resolution, and step size based on horizontal resolution.
- select_source_defaults(options, option_descriptions): Selects default values for the 'source' option based on input options.
- select_resource_defaults(options, option_descriptions): Selects default values for 'select', 'ncpus', and 'mpiprocs' options based on input options.
- find_file(pattern, path): Finds a file in the specified path that matches the given pattern.
- time_to_dhms(time_str): Converts a time string to a list of day, hour, minute, and second.
- seconds_to_dhms(seconds): Converts seconds to a list of days, hours, minutes, and seconds.
"""

import os
import fnmatch
import argparse
from datetime import datetime, timedelta
import xarray as xr
from numpy import pad

# Path to current tiegcm datafiles
TIEGCMDATA = os.environ["TIEGCMDATA"]
# Path to current tiegcm installation
TIEGCMHOME = os.environ["TIEGCMHOME"]

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
            max_ncpus = 28
            mpiprocs_default = 24
        elif o["model"] == "has":
            max_ncpus = 24
            mpiprocs_default = 24
        elif o["model"] == "ivy":
            max_ncpus = 20
            mpiprocs_default = 18
        elif o["model"] == "san":
            max_ncpus = 16
            mpiprocs_default = 12
        for on in od:
            if on == "select":
                if float(horires) == 2.5 or float(horires) == 5:
                    select = 1#4/mpiprocs_default
                if float(horires) == 1.25:
                    select = 1#8/mpiprocs_default
                if float(horires) == 0.625:
                    select = 1#12/mpiprocs_default
                select_default = int(select)
            if on == "ncpus":
                ncpus_default = max_ncpus
            if on == "mpiprocs":
                mpiprocs_default = mpiprocs_default
    return select_default,ncpus_default,mpiprocs_default

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
