#!/usr/bin/env python


"""tiegcmrun for the TIE-GCM software.

This script is performs all of the steps that are required to prepare to run a TIE-GCM simulation run. By default, this script is interactive - the user
is prompted for each decision  that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BENCH" - The user requests for a benchmark run for TIE-GCM 

"BASIC" (the default) - the user is prompted to set only a small subset of TIE-GCM
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.
"""

# Import standard modules
import argparse
import datetime
from datetime import datetime, timedelta
import json
import os
import subprocess
import fnmatch
from fractions import Fraction

# Import 3rd-party modules
#import numpy as np
from numpy import pad
import xarray as xr
from jinja2 import Template

# Import local modules
from compile import compile_tiegcm
from interpolation import interpic
from engage_solver import engage_parser, engage_run, engage_options_updater
from namelist_solver import inp_pri_date,valid_hist,inp_mxhist,inp_sechist,inp_prihist,inp_pri_out,inp_sec_out,inp_sec_date
from output_solver import create_inp_scripts, create_pbs_scripts, segment_inp_pbs
from misc import get_mtime, segment_time, valid_bench, find_file, resolution_solver, select_resource_defaults, select_source_defaults
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
        "--hidra", "-hi", action="store_true",
        help="Enable HIDRA (default: %(default)s)."
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
        if engage != None:
            od["hpc_system"]["default"] = engage["hpc_system"]
        else:
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
                        od[on]["default"] = find_file(options['simulation']['job_name']+'.nc', TIEGCMDATA)
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
        #if engage != None:
        #    print(od["SECFLDS"]["default"])
        #    od["SECFLDS"]["default"] = '["TN","UN","VN","NE","TEC","POTEN","Z","ZG","O2","O1","N2","N4S","HE","NO","NO_COOL","WN"]'
             
        temp_mode = mode
        skip_inp = []
        segment = None
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
                PRIHIST = inp_prihist(PRISTART,PRISTOP, segment)
                od["PRIHIST"]["default"] =  PRIHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                PRIHIST = [int(i) for i in o[on].split()]
                MXHIST_PRIM_set ,MXHIST_PRIM_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], PRIHIST, od["MXHIST_PRIM"]["warning"], segment)
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
                SECHIST = inp_sechist(PRISTART,PRISTOP, segment)
                od["SECHIST"]["default"] =  SECHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters)
                SECHIST = [int(i) for i in o[on].split()]
                MXHIST_SECH_set ,MXHIST_SECH_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], SECHIST, od["MXHIST_SECH"]["warning"],segment)
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
            if engage != None:
                od['group_list']['default'] = engage['group_list']
                od['queue']['default'] = engage['queue']
                od['walltime']['default'] = engage['walltime']
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
    hidra = args.hidra
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
    submit_all_jobs_script = ''
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
        if engage != None:
            with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
                option_descriptions = json.load(f)
            options = engage_options_updater(options, args.engage, option_descriptions)
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
    json_path = f"{workdir}/{run_name}.json"
    
    tiegcmdata = options["model"]["data"]["tgcmdata"]
    horires = options["model"]["specification"]["horires"]
    vertres = options["model"]["specification"]["vertres"]
    zitop = options["model"]["specification"]["zitop"]
    
    if options.get("inp") == None:
        input_file_generatred = True
    else:
        input_file_generatred = False
    if args.onlycompile == True:
        compile_tiegcm(options, debug, coupling, hidra)
    elif args.engage != None:
        options_coupling,standalone_pbs_files,coupling_inp_files = engage_run(options, debug, coupling, args.engage)
        return (options_coupling,standalone_pbs_files,coupling_inp_files)
    else:
        if args.compile == True:
            compile_tiegcm(options, debug, coupling, hidra)
        if options["model"]["specification"]["segmentation"] == False:
            if options["model"]["data"]["input_file"] == None or not os.path.isfile(options["model"]["data"]["input_file"]):
                if input_file_generatred == False:
                    if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
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

            if os.path.exists(json_path):
                if not clobber:
                    raise FileExistsError(f"Options file {json_path} exists!")
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(options, f, indent=JSON_INDENT)

            if options["simulation"]["hpc_system"] != "linux":
                pbs_script = create_pbs_scripts(options,run_name,None)
        else:
            if not os.path.isfile(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'):
                in_prim = options["inp"]["SOURCE"]
                out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                options["inp"]["SOURCE"] = out_prim
                interpic (in_prim,float(horires),float(vertres),float(zitop),out_prim)
            else:
                out_prim = f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc'
                options["inp"]["SOURCE"] = out_prim
                print(f'{options["model"]["data"]["workdir"]}/{run_name}_prim.nc exists')
            inp_files, pbs_files,log_files, pristart_times, pristop_times = segment_inp_pbs(options, run_name, pbs=True)
            init_inp = inp_files[0]    
            init_pbs = pbs_files[0]
            options["model"]["data"]["input_file"] = init_inp

            # Create a single script which will submit all of the PBS jobs in order.
            os.chdir(workdir)
            submit_all_jobs_script = (f"{options['simulation']['job_name']}_pbs.sh")
            with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
                cmd = f"#!/bin/bash\n"
                f.write(cmd)
                cmd = f"# TIEGCM Jobs\n"
                f.write(cmd)
                tiegcm_pbs = pbs_files[0]
                cmd = f"tiegcm_job_id=`qsub {tiegcm_pbs}`\n"
                f.write(cmd)
                cmd = "echo $tiegcm_job_id\n"
                f.write(cmd)
                for tiegcm_pbs in pbs_files[1:]:
                    cmd = "old_tiegcm_job_id=$tiegcm_job_id\n"
                    f.write(cmd)
                    cmd = f"tiegcm_job_id=`qsub -W depend=afterok:$old_tiegcm_job_id {tiegcm_pbs}`\n"
                    f.write(cmd)
                    cmd = "echo $tiegcm_job_id\n"
                    f.write(cmd)
            os.chmod(submit_all_jobs_script, 0o755)

    if args.execute == True and options["simulation"]["hpc_system"] != "linux":
        if args.compile == False:
            if find_file(options["model"]["data"]["modelexe"], execdir) == None and os.path.exists(options["model"]["data"]["modelexe"]) == False :
                print(f'{RED}Unable to find executable in {execdir}{RESET}')
                exit(1)
        try:
            if submit_all_jobs_script == '':
                result = subprocess.run(['qsub', pbs_script], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                job_id = result.stdout.strip()
                print(f'Job submitted successfully. Job ID: {job_id}')
            else:
                result = subprocess.run(['./'+submit_all_jobs_script], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                print(f'Jobs submitted successfully')
        except subprocess.CalledProcessError as e:
            print(f'{YELLOW}Error submitting job: {e.stderr}{RESET}')
            print(f"{YELLOW}Check PBS script for errors{RESET}")
            if submit_all_jobs_script == '':
                print(f"To submit job use command {YELLOW}qsub {pbs_script}{RESET}")
            else:
                print(f"To submit job use command {YELLOW} ./{submit_all_jobs_script}{RESET}")
    
    elif args.onlycompile == False and options["simulation"]["hpc_system"] != "linux":
        print(f"{YELLOW}Execute is set to false{RESET}")
        if submit_all_jobs_script == '':
                print(f"To submit job use command {YELLOW}qsub {pbs_script}{RESET}")
        else:
            print(f"To submit job use command {YELLOW} ./{submit_all_jobs_script}{RESET}")
    elif args.onlycompile == False:
        print(f"{YELLOW}HPC System is set to linux{RESET}")
        print(f"{YELLOW}To run the model use command{RESET} mpirun {options['model']['data']['modelexe']} {options['model']['data']['input_file']}")

if __name__ == "__main__":
    """Begin tiegcmrun program."""
    tiegcmrun()
