"""
namelist_solver.py for the TIEGCMrun software.

This script helps in solving for the correct input parameter for the namelist files for running the TIEGCM model. 

Functions included:

- inp_pri_date(start_date_str, stop_date_str): Converts start and stop date strings to datetime objects and extracts relevant information.
- valid_hist(start_time, stop_time): Calculates valid divisions for a given time range.
- inp_mxhist(start_time, stop_time, x_hist, mxhist_warn, segment=None): Calculates the MXHIST value based on the given start and stop times, and x_hist values.
- inp_sechist(SECSTART, SECSTOP, segment=None): Determines the value of SECHIST based on the given SECSTART and SECSTOP.
- inp_prihist(PRISTART, PRISTOP, segment=None): Calculates the PRIHIST list based on the PRISTART and PRISTOP values.
- inp_pri_out(start_time, stop_time, PRIHIST, MXHIST_PRIM, pri_files, histdir, run_name): Generates the output file names and the total number of files based on the given parameters.
- inp_sec_out(start_time, stop_time, SECHIST, MXHIST_SECH, sec_files, histdir, run_name): Calculates the output file names and the total number of files for the secondary history output.
- inp_sec_date(start_time, stop_time, SECHIST): Calculates the SECSTART and SECSTOP values based on the given start_time, stop_time, and SECHIST.
"""


from datetime import datetime, timedelta
from math import ceil


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

def inp_mxhist(start_time, stop_time, x_hist, mxhist_warn, segment = None):
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
    if segment == None:
        if mxhist_day >= 1:
            mxhist_warn = (mxhist_warn + "\n" if mxhist_warn is not None else "") + f"For a Daily output set MXHIST to {int(mxhist_day)}"
        if mxhist_hour >= 1:
            mxhist_warn = mxhist_warn +  f"\nFor a Hourly output set MXHIST to {int(mxhist_hour)}"
        if mxhist_min >= 1:
            mxhist_warn = mxhist_warn + f"\nFor a Minutely output set MXHIST to {int(mxhist_min)}"
        MXHIST = mxhist_day
    else:
        segment_seconds = (segment[0] * 86400) + (segment[1] * 3600) + (segment[2] * 60) + segment[3]
        MXHIST = segment_seconds/step_seconds
        mxhist_warn = f"MXHIST minimum = 1, maximum = {int(MXHIST)} for segment run."
    return(int(MXHIST), mxhist_warn)

def inp_sechist(SECSTART, SECSTOP, segment = None):
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
    
    if segment != None:
        SECHIST = [1 if x != 0 else 0 for x in segment]

    return SECHIST

def inp_prihist(PRISTART, PRISTOP, segment = None):
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

    if segment != None:
        PRIHIST = [1 if x != 0 else 0 for x in segment]
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
                OUTPUT = OUTPUT = f"'{histdir}/{run_name}_temp_{'{:02d}'.format(pri_files)}.nc' , '{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files+1)}.nc'"
        else:
            PRIM_0 = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
            PRIM_N = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files_n)}.nc"
            OUTPUT = f"'{PRIM_0}','to','{PRIM_N}','by','1'"
    else:
        if number_of_files == 1:
            OUTPUT = OUTPUT = f"'{histdir}/{run_name}_temp_{'{:02d}'.format(pri_files)}.nc' , '{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files+1)}.nc'"
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
    sec_files_start = sec_files + 1  # Start numbering from next file
    sec_files_end = sec_files_start + number_of_files  # End numbering based on number of files
    if number_of_files == 1 or number_of_files == 0:
        # If only one file is being generated
        SECOUT = f"'{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_start)}.nc'"
        sec_files_end = sec_files_start
    else:
        # If multiple files are being generated
        SECH_0 = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_start)}.nc"
        SECH_N = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_end)}.nc"
        SECOUT = f"'{SECH_0}','to','{SECH_N}','by','1'"

    # Return the new sec_files value for subsequent calls
    return SECOUT, sec_files_end

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
