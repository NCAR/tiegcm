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

from misc import seconds_to_dhms, resolution_solver, find_file, select_resource_defaults, get_mtime, select_source_defaults
from output_solver import segment_inp_pbs
from interpolation import interpic
from namelist_solver import inp_pri_date


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
    elif hpc_system == 'pleiades':
        group_list = o['group_list']
    walltime = o['walltime']
    modules = o['modules']

    o = engage_parameters["coupling"]
    
    gr_warm_up_time = int(o['gr_warm_up_time'])
    gcm_spin_up_time = int(o['gcm_spin_up_time'])
    conda_env = o['conda_env']
    
    start_date = get_engage_start_time(coupled_start_date,gr_warm_up_time+gcm_spin_up_time)
    
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
    if hpc_system == 'pleiades':
        eo['group_list'] = group_list
    else:
        eo['group_list'] = None
    eo['queue'] = queue
    
    eo['walltime'] = walltime
    eo['modules'] = modules
    eo['conda_env'] = conda_env
    
    if hpc_system == "derecho":
        eo['job_priority'] = job_priority
    elif hpc_system == 'pleiades':
        eo["model"] = "bro"
    
    eo['skip']= ['group_list','job_name','hpc_system','horires','parentdir','vertres', 'mres', 'input_file', 'LABEL','start_time','stop_time','secondary_start_time','secondary_stop_time','segment' ,'SOURCE_START','PRIHIST','MXHIST_PRIM','SECHIST','MXHIST_SECH','account_name','project_code','queue','job_priority','model','walltime']
    
    return engage_options

def engage_options_updater(options, engage_options, option_descriptions):
    # General options for the simulation
    o = options["simulation"] 
    o["job_name"] = engage_options["job_name"]
    run_name = o["job_name"]
    o["hpc_system"] = engage_options["hpc_system"]
    # Data Options
     
    o = options["model"]["data"]
    o["parentdir"] = engage_options["parentdir"]
    o["execdir"] = o["parentdir"]
    o["workdir"] = o["parentdir"]
    o["histdir"] = o["parentdir"]
    # Specification Options
    o = options["model"]["specification"]
    o["horires"] = engage_options["horires"]
    horires = o["horires"]
    vertres, mres, nres_grid, STEP = resolution_solver(o["horires"])
    if o.get("vertres") is None:
        o["vertres"] = vertres
    if o.get("mres") is None:
        o["mres"] = mres
    if o.get("nres_grid") is None:
        o["nres_grid"] = nres_grid
    zitop = o["zitop"]
    # INP options
    o = options["inp"]
    if o.get("STEP") is None:
        o["STEP"] = STEP
    o["start_time"] = engage_options["start_time"]
    o["stop_time"] = engage_options["stop_time"]
    o["segment"] = " ".join(map(str, engage_options["segment"]))
    options_temp = copy.deepcopy(options)
    if o.get("SOURCE") is None:
        print("No SOURCE file specified, creating a new one.")
        o["SOURCE"] = select_source_defaults(options_temp, option_descriptions)  
        """
        if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/tiegcm_standalone/{run_name}-tiegcm-standalone_temp.nc'):
            in_prim = source
            out_prim = f'{options["model"]["data"]["workdir"]}/tiegcm_standalone/{run_name}-tiegcm-standalone_temp.nc'
            o["SOURCE"] = out_prim
            interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
        """
    if o.get("SOURCE_START") is None:
        o["SOURCE_START"] =  " ".join(map(str, get_mtime(options["inp"]["SOURCE"])[0]))
    START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(o["start_time"], o["stop_time"])
    if o.get("START_YEAR") is None:
        o["START_YEAR"] = START_YEAR
    if o.get("START_DAY") is None:
        o["START_DAY"] = START_DAY
    # PBS options
    o = options["job"]
    o["account_name"] = engage_options["account_name"]
    hpc_platform = options["simulation"]["hpc_system"]
    if hpc_platform == "derecho":
        if o.get("mpi_command") is None:
            o["mpi_command"] = "mpirun"
        o['queue'] = engage_options['queue']
        o['job_priority'] = engage_options['job_priority']
        o['walltime'] = engage_options['walltime']
    elif hpc_platform == "pleiades":
        o['group_list'] = engage_options['group_list']
        if o.get("mpi_command") is None:
            o["mpi_command"] = "mpiexec_mpt"
        o['queue'] = engage_options['queue']
        o['walltime'] = engage_options['walltime']
    on = options["job"]["resource"] = {}
    if hpc_platform == "derecho":
        select_default,ncpus_default,mpiprocs_default = select_resource_defaults(options,option_descriptions)
        if on.get("select") is None:
            on["select"] = select_default
        if on.get("ncpus") is None:
            on["ncpus"] = ncpus_default
        if on.get("mpiprocs") is None:
            on["mpiprocs"] = mpiprocs_default
    elif hpc_platform == "pleiades":
        if on.get("model") is None:
            on["model"] = "bro"
        select_default,ncpus_default,mpiprocs_default = select_resource_defaults(options,option_descriptions)
        if on.get("select") is None:
            on["select"] = select_default
        if on.get("ncpus") is None:
            on["ncpus"] = ncpus_default
        if on.get("mpiprocs") is None:
            on["mpiprocs"] = mpiprocs_default
        if on.get("moduledir") is None:
            options["job"]["moduledir"] = option_descriptions["job"]["pleiades"]["moduledir"]["default"]
        if on.get("local_modules") is None:
            options["job"]["local_modules"] = option_descriptions["job"]["pleiades"]["local_modules"]["default"]
        if on.get("other_job") is None:
            options["job"]["other_job"] = option_descriptions["job"]["pleiades"]["other_job"]["default"]
    
    o["nprocs"] = int(on["select"]) * int(on["mpiprocs"])
    o["modules"] = engage_options["modules"]
    o["project_code"] = engage_options["project_code"]
    
    return options


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
    horires_standalone = float(engage["horires"])
    options_standalone["inp"]["stop_time"] = engage["coupled_start_time"]
    options_standalone["inp"]["PRIHIST"] = '1 0 0 0'
    options_standalone["inp"]["MXHIST_PRIM"] = 1
    options_standalone["inp"]["SECHIST"] = '0 1 0 0'
    options_standalone["inp"]["MXHIST_SECH"] = 24
    options_standalone["inp"]["segment"] = '1 0 0 0'
    options_standalone["inp"]["OPDIFFCAP"] = '2e9'
    options_standalone["inp"]["OPDIFFRATE"] = '0.3'
    options_standalone["inp"]["OPDIFFLEV"] = '7'
    options_standalone["inp"]["OPFLOOR"] = '3000'
    options_standalone["inp"]["OPRATE"] = '0.3'
    options_standalone["inp"]["OPLEV"] = '7'
    options_standalone["inp"]["OPLATWIDTH"] = '20'
    options_standalone["inp"]["TE_CAP"] = '8000'
    options_standalone["inp"]["TI_CAP"] = '8000'
    options_standalone["inp"]["GSWM_MI_DI_NCFILE"] = find_file(f'*gswm_diurn_{horires_standalone}d_99km*', TIEGCMDATA)
    options_standalone["inp"]["GSWM_MI_SDI_NCFILE"] = find_file(f'*gswm_semi_{horires_standalone}d_99km*', TIEGCMDATA)
    options_standalone["inp"]["GSWM_NM_DI_NCFILE"] = find_file(f'*gswm_nonmig_diurn_{horires_standalone}d_99km*', TIEGCMDATA)
    options_standalone["inp"]["GSWM_NM_SDI_NCFILE"] = find_file(f'*gswm_nonmig_semi_{horires_standalone}d_99km*', TIEGCMDATA)
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
    options_coupling["inp"]["OPDIFFCAP"] = '2e9'
    options_coupling["inp"]["OPDIFFRATE"] = '0.3'
    options_coupling["inp"]["OPDIFFLEV"] = '7'
    options_coupling["inp"]["OPFLOOR"] = '3000'
    options_coupling["inp"]["OPRATE"] = '0.3'
    options_coupling["inp"]["OPLEV"] = '7'
    options_coupling["inp"]["OPLATWIDTH"] = '20'
    options_coupling["inp"]["TE_CAP"] = '8000'
    options_coupling["inp"]["TI_CAP"] = '8000'
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
