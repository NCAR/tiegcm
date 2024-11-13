"""
engage_solver.py for the TIEGCMrun software.

This script helps in solving for the correct input and pbs parameters for running the TIEGCM model in coupled mode with Engage. 

Functions included:

- gamres_to_res(gamres): Converts the GAMERA grid type to horizontal resolution values.
- engage_parser(engage_parameters): Parses the engage.json file and returns the options dictionary.
- get_engage_start_time(datetime_str, seconds): Calculates the start time for the Engage run by subtracting the spin-up time from the coupled start date.
- engage_run(options, debug, coupling, engage): Prepares and runs the TIEGCM model in both standalone and coupled modes, generating the necessary input and PBS files.
"""

import os
import json
import copy
from datetime import datetime, timedelta

from misc import seconds_to_dhms, resolution_solver, find_file, select_resource_defaults
from output_solver import segment_inp_pbs
from interpolation import interpic


# Path to current tiegcm datafiles
TIEGCMDATA = os.environ["TIEGCMDATA"]
# Path to current tiegcm installation
TIEGCMHOME = os.environ["TIEGCMHOME"]
# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(TIEGCMHOME, "tiegcmrun")
OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")

def gamres_to_res(gamres):
    "D", "Q", "O", "H"
    if gamres == "D":
        return 2.5 , 2.5
    elif gamres == "Q":
        return 2.5 , 1.25
    elif gamres == "O":
        return 1.25 , 0.625
    elif gamres == "H":
        return 1.25 ,0.625

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
    if hpc_system == "derecho":
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
    
    eo['walltime'] = walltime
    eo['modules'] = modules
    eo['conda_env'] = conda_env
    
    if hpc_system == "derecho":
        eo['job_priority'] = job_priority
    elif hpc_system == 'pleiades':
        eo["model"] = "bro"
    
    eo['skip']= ['job_name','hpc_system','horires','parentdir','vertres', 'mres', 'input_file', 'LABEL','start_time','stop_time','secondary_start_time','secondary_stop_time','segment' ,'SOURCE_START','PRIHIST','MXHIST_PRIM','SECHIST','MXHIST_SECH','account_name','project_code','queue','job_priority','model','walltime']
    
    return engage_options

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
