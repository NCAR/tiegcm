"""
output_solver.py for the TIEGCMrun software.

This script generates input and PBS scripts for running the TIEGCM model. It handles the creation of segmented input and PBS scripts based on the provided model options and time segments.

Functions included:

- create_pbs_scripts(options, run_name, segment_number): Creates PBS scripts for running the TIEGCM model.
- create_inp_scripts(options, run_name, segment_number): Creates input scripts for running the TIEGCM model.
- segment_inp_pbs(options, run_name, pbs, engage_options=None): Segments the input and PBS scripts based on the provided options and time segments.
"""


import os
import json
import copy
from jinja2 import Template

from misc import segment_time, resolution_solver
from namelist_solver import inp_pri_date, inp_pri_out, inp_sec_date, inp_sec_out


JSON_INDENT = 4
# Path to current tiegcm datafiles
TIEGCMDATA = os.environ["TIEGCMDATA"]
# Path to current tiegcm installation
TIEGCMHOME = os.environ["TIEGCMHOME"]
# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(TIEGCMHOME, "tiegcmrun")
OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")
# Path to template .inp file.
INP_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.inp")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

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
                    segment_PRIHIST = [int(i) for i in PRIHIST.split()]
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
                if segment_number == last_segment_time and engage_options != None:
                    interpolation_script = os.path.join(segment_options["model"]["data"]["workdir"],f'tiegcm_resolution_upscale.py')
                    with open(interpolation_script, "w", encoding="utf-8") as f:
                        f.write("import sys\n")
                        f.write(f"sys.path.append('{TIEGCMHOME}/tiegcmrun')\n")
                        f.write("import interpolation\n")
                        horires_coupled = engage_options["horires_coupled"]
                        vertres_coupled, mres_coupled, nres_grid_coupled, STEP_coupled = resolution_solver(horires_coupled,engage_options)
                        SOURCE_coupling = os.path.join(os.path.dirname(segment_options["model"]["data"]["workdir"]),f'{engage_options["job_name"]}_prim.nc')
                        input_standalone = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
                        f.write(f"interpolation.interpic('{input_standalone}',{float(horires_coupled)},{float(vertres_coupled)},{float(segment_options['model']['specification']['zitop'])},'{SOURCE_coupling}')\n")
                        if options["simulation"]["hpc_system"] == "derecho":
                            interpolation_pbs = [f'conda activate {engage_options["conda_env"]}',f'python {interpolation_script}']
                        elif options["simulation"]["hpc_system"] == "pleiades":
                            interpolation_pbs =  [f'source activate {engage_options["conda_env"]}',f'python {interpolation_script}']
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
