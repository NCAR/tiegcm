#!/usr/bin/env python


"""tiegcmrun: interactively configure, build and submit a TIE-GCM run.

--mode (BENCH/BASIC/INTERMEDIATE/EXPERT) sets which options are prompted; the rest take defaults.
"""

import argparse
import copy
import json
import os
import re
import socket
import subprocess
import sys
from fractions import Fraction

from compile import compile_tiegcm
from interpolation import interpic
from engage_solver import engage_parser, engage_run, engage_options_updater, apply_coupled_defaults, COUPLED_DEFAULT_KEYS, COUPLED_HISTORY_KEYS
from namelist_solver import inp_pri_date,inp_mxhist,inp_sechist,inp_prihist,inp_pri_out,inp_sec_out,inp_sec_date,parse_segment,forcing_missing,parse_run_datetime,cadence_step,hms_str,hms_differs,source_start_hms_error,source_time_of_day_error,resolve_source_start,forcing_refusal,date_order_problem,year_boundary,lbc_other_set
from output_solver import create_inp_scripts, create_pbs_scripts, segment_jobs, write_segment_jobs, run_file, render_pbs, CONFIG_DIR, apply_system_config, linux_run_commands, queue_walltime_default
from jobgen import JobGen, spell, cdpath_reset
from misc import get_mtime, segment_time, valid_bench, find_file, resolution_solver, mres_to_nres_grid, default_make_file, select_resource_defaults, select_source_defaults, select_latest_gpi, file_covers, resolve_benchmark_file, generate_gpi_file, generate_imf_file, as_list, require_env, tiegcm_env, in_execdir, GSWM_PATTERNS, gswm_allowed, he_coefs_pattern, unknown_other_input_keys, other_input_error, secflds_list, ascii_quotes
from replay_solver import (prepare_replay, parse_rederive, rederive_pop, rederive_groups_text,
                           old_format_problem, validate_options, read_options_json, section_problem,
                           options_key_schema, unknown_keys_problem)
import misc
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
RESET = '\033[0m'


DESCRIPTION = "Interactive script to prepare a TIEGCM model run."

# engage refuses a tiegcmrun whose major version differs; bump the major on an incompatible change.
API_VERSION = "2.1.0"

JSON_INDENT = 4

PARAMETERS_JSON = "tiegcmrun_parameters.json"

SUPPORT_FILES_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")

BENCHMARKS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, 'benchmarks.json')


MODES = ["BENCH", "BASIC", "INTERMEDIATE", "EXPERT"]


def _require_tiegcm_env():
    require_env("TIEGCMHOME")
    require_env("TIEGCMDATA")


def _import_own_tui_form():
    """Load this package's tui_form.py by path; under engage a bare import could bind makeitso's."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "tiegcm_tui_form", os.path.join(here, "tui_form.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def use_tui(no_tui):
    """True unless --no-tui or stdin/stdout is not a terminal."""
    return not no_tui and sys.stdin.isatty() and sys.stdout.isatty()


def create_command_line_parser():
    """Return the tiegcmrun command-line argument parser."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--clobber", action="store_true",
        help="Overwrite existing parameters JSON and .inp/.pbs files (default: False)."
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
        "--debug-build", action="store_true",
        help="With --compile/--onlycompile: build a debug executable (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default=None, type=str.upper, choices=MODES, metavar="MODE",
        help="BENCH|BASIC|INTERMEDIATE|EXPERT (default: BASIC, BENCH with -bench)."
    )
    parser.add_argument(
        "--no-tui", action="store_true",
        help="Ask the run options as linear prompts instead of the wizard."
    )
    parser.add_argument(
        "--tui", action="store_true",
        help="The wizard (the default; kept for old command lines)."
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
        help="Path to JSON file of options (default: %(default)s)."
    )
    parser.add_argument(
        "--rederive", default=None, metavar="GROUP[,GROUP]",
        help="With -o: re-derive these groups; '--rederive list' names them."
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
        help="Benchmark run name (default: %(default)s)."
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Validate the options and stop before anything is written."
    )
    parser.add_argument(
        "--engage", default=None,
        help="Inline JSON string of engage options (default: %(default)s)."
    )
    return parser

def get_run_option(name, description, mode="BASIC", skip_parameters=[], run_dir=None):
    """Prompt for one run option and return its value; options above `mode` return the default.

    A typed relative path resolves against run_dir (None: as typed).
    """
    level = description["LEVEL"]
    prompt = description.get("prompt", "")
    default = description.get("default", None)
    valids = description.get("valids", None)
    var_description = description.get("description", None)
    warning = description.get("warning", None)
    fourvar_variables = misc.keys_of_type("dhms")

    if mode == "BENCH" and level in ["BASIC","INTERMEDIATE", "EXPERT"]:
        if name in fourvar_variables and default is not None:
            return  ' '.join(map(str, default))
        else:
            return default
    # Never prompted in any mode: a different answer would desync the coupled deck.
    if name in skip_parameters:
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
    file_variables = misc.INP_FILE_KEYS
    if valids is not None: 
        if name == "vertres":
            vs = "|".join(map(lambda x: str(Fraction(x)), valids))
            prompt += f" ({vs})"            
        else:
            vs = "|".join(map(str,valids))
            prompt += f" ({vs})"
    elif description.get("choices"):
        # suggestions only: any other value is accepted
        prompt += f" (your projects: {'|'.join(map(str, description['choices']))})"

    if default is not None:
        if name == "vertres":
            prompt += f" [{GREEN}{str(Fraction(default))}{RESET}]"
        else:
            prompt += f" [{GREEN}{default}{RESET}]"

    ok = False
    option_value = ""
    secflds_typed = None
    while not ok:
        if name in ("other_input", "other_pbs","other_job","job_chain"):
            prompt = og_prompt
            prompt += f" [{GREEN}{default}{RESET}]"
            temp_value = input(f"{prompt} / ENTER to go next: ")
            if temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                option_value = default
                ok = True
            else:
                # A key that is not a &tgcm_input member would stop the model's namelist read.
                entries = [x for x in as_list(temp_value) if x is not None]
                unknown = unknown_other_input_keys(entries) if name == "other_input" else []
                if unknown:
                    print(f'{YELLOW}{other_input_error(unknown)}{RESET}')
                    continue
                default = list(default or []) + entries
                option_value = default
        elif name == "SECFLDS":
            prompt = og_prompt
            prompt += f" [{GREEN}{default}{RESET}]"
            temp_value = input(f"{prompt} / ENTER to go next: ")
            if temp_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
            elif temp_value == 'none' or temp_value == 'None':
                option_value = json.loads('[null]')
                ok = True
            elif temp_value == "":
                option_value = default
                ok = True
            else:
                # Typed names replace the default list: the model stops on duplicate SECFLDS names.
                typed = (secflds_typed or []) + [temp_value]
                try:
                    names, notes = secflds_list(
                        [n for line in typed for n in as_list(ascii_quotes(line), tokens=True)])
                except ValueError as e:
                    print(f'{RED}{e}{RESET}')
                    print(f'{YELLOW}Re-enter the field names (ENTER keeps {default}).{RESET}')
                    continue
                if notes:
                    print(f'{YELLOW}SECFLDS: {"; ".join(dict.fromkeys(notes))}{RESET}')
                secflds_typed = typed
                default = names
                option_value = default
        elif name == "input_file":
            # The user's own deck: no TIEGCMDATA lookup.
            option_value = input(f"{prompt}: ")
            if option_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
                continue
            if option_value == "":
                option_value = default
            elif option_value in ("none", "None"):
                option_value = None
            else:
                path = os.path.abspath(misc.anchor_path(option_value, run_dir) if run_dir
                                       else os.path.expanduser(option_value))
                if not os.path.isfile(path):
                    print(f'{YELLOW} Custom input file {option_value} does not exist '
                          f'({"a directory" if os.path.isdir(path) else "no such file"}). '
                          f'Give the path to an existing .inp, or ENTER to build one.{RESET}')
                    continue
                option_value = path
            ok = True
        elif name in file_variables:
            option_value = input(f"{prompt}: ")
            if option_value == "":
                option_value = default
                if name == "SOURCE" and option_value == None:
                    continue
                else:
                    ok = True
            elif option_value == 'none' or option_value == 'None':
                option_value = None
                ok = True
            elif option_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
                continue
            elif option_value.strip().lower() == "gen" and name in ("GPI_NCFILE", "IMF_NCFILE"):
                option_value = "gen"  # handled by the caller
                ok = True
            elif "$" in option_value:
                # resolved by the job's shell: no existence check
                option_value = option_value.strip()
                ok = True
            elif option_value != None:
                # relative to the run folder; a bare name not found there is looked up under TIEGCMDATA
                path = misc.anchor_path(option_value, run_dir) if run_dir else option_value
                if os.path.isfile(path) == False:
                    if os.path.isdir(path) == True:
                        print(f'{YELLOW} {option_value} is a directory. Please provide a file path.{RESET}')
                        continue
                    else:
                        file_path = find_file(option_value, tiegcm_env("TIEGCMDATA"))
                        if file_path == None:
                            print(f'{YELLOW} Unable to find {option_value} in {tiegcm_env("TIEGCMDATA")}.\n Give path to file as an alternative.{RESET}')
                            continue
                        else:
                            print(f'File Found: {file_path}')
                            option_value = str(file_path)
                            ok = True
                else:
                    option_value = str(path)
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
                else:
                    temp_value = temp_value.replace("'", "")
                    temp_array.extend(s.replace(" ", "") for s in temp_value.split()) 
                try:
                    temp_array = [int(i) for i in temp_array]
                except ValueError:
                    print(f'{YELLOW}Invalid value for option {name}: {temp_value!r} '
                          f'(give 4 integers: Day Hour Min Sec){RESET}')
                    continue
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
                # only the segment length is optional ('none' = no segmentation)
                if name != "segment":
                    print(f'{YELLOW}{name} is required: give 4 integers (Day Hour Min Sec), or ENTER '
                          f'for the default.{RESET}')
                    continue
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
            option_value = input(f"{prompt}: ")

            # A list-valued option must stay a list: the .pbs renderer would iterate a str by character.
            if isinstance(default, list):
                if option_value == "?":
                    print(f'{YELLOW}{var_description}{RESET}')
                    continue
                option_value = default if option_value == "" else as_list(option_value)
                ok = True
                continue

            # Only a typed relative path is resolved against the run folder, not a default.
            typed = option_value not in ("", "none", "None")
            if option_value == "":
                option_value = default
            elif option_value == 'none' or option_value == 'None':
                option_value = None
            elif option_value == "?":
                print(f'{YELLOW}{var_description}{RESET}')
                continue
            try:
                if name == "vertres":
                    if valids is not None and float(Fraction(option_value)) not in valids:
                        print(f"Invalid value for option {name}: {option_value}!")
                        continue
                    else:
                        option_value = str(float(Fraction(option_value)))
                elif name in _REQUIRED_DATES and option_value is None:
                    parse_run_datetime(option_value, name)     # raises: 'none' is refused
                elif name in _REQUIRED and option_value is None:
                    raise ValueError(f"{name} is required")
                elif option_value is not None:
                    option_value = misc.coerce(name, option_value, description,
                                               run_dir if typed else None)
            except misc.ValidsError:
                print(f"Invalid value for option {name}: {option_value}!")
                continue
            except (ValueError, TypeError, ZeroDivisionError) as e:
                print(f"{YELLOW}{e}{RESET}")
                continue
            ok = True
    return option_value

# Answers that cannot be 'none'.
_REQUIRED_DATES = ("start_time", "stop_time", "secondary_start_time", "secondary_stop_time")
_REQUIRED = ("horires", "zitop", "mres", "ELECTRON_HEATING", "MXHIST_PRIM", "MXHIST_SECH")


# Keys derived from the run window; a segmented run derives them per segment and ignores an edit.
DERIVED_INP_KEYS = ("START_YEAR", "START_DAY", "PRISTART", "PRISTOP", "SECSTART", "SECSTOP",
                    "OUTPUT", "SECOUT")


def _derived_note(name, value, default, segmented=False):
    """Note an answer that differs from the derived default; True means keep the default."""
    shown = " ".join(map(str, default)) if isinstance(default, (list, tuple)) else default
    if value is None or shown is None or str(value) == str(shown):
        return False
    if segmented:
        print(f"{YELLOW}NOTE: {name} {value} not used: a segmented run derives it per segment.{RESET}")
        return True
    print(f"{YELLOW}NOTE: {name} {value} differs from the derived {shown}; written as given.{RESET}")
    return False


def _prompted(name, description, mode, skip_parameters):
    """True when get_run_option(name, ...) asks the user (False: it returns the default)."""
    level = description["LEVEL"]
    if name in skip_parameters or (mode == "BENCH" and level != "BENCH"):
        return False
    if mode == "BASIC" and level in ("INTERMEDIATE", "EXPERT"):
        return False
    if mode == "INTERMEDIATE" and level == "EXPERT":
        return False
    return True


def _reprompt_while(name, description, mode, skip_parameters, value, problem):
    """Re-ask `name` while problem(value) returns an error; raises if `name` is not prompted."""
    err = problem(value)
    while err:
        if not _prompted(name, description, mode, skip_parameters):
            raise ValueError(f"{RED}{err}{RESET}")
        print(f"{RED}{err}{RESET}")
        value = get_run_option(name, description, mode, skip_parameters)
        err = problem(value)
    return value


_date_order_problem = date_order_problem


_FORCING_REPROMPTS = 3

def _complete_forcing(o, od, mode, skip_parameters):
    """Re-prompt the solar/high-latitude forcing fields input.F requires; raise if still unset.

    KP is offered first when POWER/CTPOTEN are missing (it stands in for both).
    """
    for _ in range(_FORCING_REPROMPTS):
        missing = forcing_missing(o)
        if not missing:
            return
        names = [f for f, _ in missing]
        blocked = [f for f in names if f in skip_parameters]
        if blocked:
            break
        print(f"{RED}Incomplete solar/high-latitude forcing (the model would stop):{RESET}")
        for _, why in missing:
            print(f"{RED}  - {why}{RESET}")
        try:
            if ("POWER" in names or "CTPOTEN" in names) and "KP" not in skip_parameters:
                o["KP"] = get_run_option("KP", od["KP"], mode, skip_parameters)
                names = [f for f, _ in forcing_missing(o)]
            for fld in names:
                o[fld] = get_run_option(fld, od[fld], mode, skip_parameters)
        except EOFError:
            break
    missing = forcing_missing(o)
    if missing:
        raise ValueError(f"{RED}{forcing_refusal(missing)}{RESET}")

_hms_str = hms_str
_hms_differs = hms_differs

def _source_start_candidates(o, od, source_mtimes, pristart, skip_parameters, run_dir=None):
    """Return the SOURCE histories at PRISTART's time of day (input.F requires equal hr/min/sec).

    With none, re-prompts SOURCE (updating o['SOURCE']); keeping the same file raises.
    """
    matching = [m for m in source_mtimes if list(m[1:4]) == list(pristart[1:4])]
    while not matching:
        message, refusal = source_time_of_day_error(o["SOURCE"], source_mtimes, pristart, o["start_time"])
        print(f"{RED}{message}{RESET}")
        previous = o["SOURCE"]
        try:
            o["SOURCE"] = get_run_option("SOURCE", od["SOURCE"], "BASIC", skip_parameters, run_dir=run_dir)
        except EOFError:
            raise ValueError(refusal) from None
        if o["SOURCE"] in (None, previous):
            raise ValueError(refusal)
        source_mtimes = get_mtime(o["SOURCE"])
        matching = [m for m in source_mtimes if list(m[1:4]) == list(pristart[1:4])]
    return matching

def prompt_user_for_run_options(args):
    """Prompt for the run options, group by group, and return the options dict."""

    input_build_skip = False
    pbs_build_skip = False
    base_skip = False
    skip_parameters = []
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
        # fields that do not affect the build
        skip_parameters = ['input_file','log_file','job_name','modeldir','parentdir','tgcmdata','segmentation']
        if coupling == True:
            skip_parameters.append('modelexe')
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = apply_system_config(json.load(f), socket.gethostname())
    with open(BENCHMARKS_FILE, "r", encoding="utf-8") as f:
        benchmarks_options = json.load(f)
    options = {}

    if base_skip == False:
        o = options["simulation"] = {}
        od = option_descriptions["simulation"]
        if benchmark != None:
            od["job_name"]["default"] = benchmark
        elif engage != None:
            od["job_name"]["default"] = engage["job_name"]
        if engage != None:
            od["hpc_system"]["default"] = engage["hpc_system"]
        for on in ["job_name", "hpc_system"]:
            o[on] = get_run_option(on, od[on], mode, skip_parameters)

    options["model"] = {}

    options["model"]["data"] = {}
    o = options["model"]["data"]
    temp_mode = mode
    od = option_descriptions["model"]["data"]
    od["modeldir"]["default"] = tiegcm_env("TIEGCMHOME")
    # asked before the run folder is known: a relative answer is resolved against it below
    o["modeldir"] = get_run_option("modeldir", {**od["modeldir"], "type": "path"}, mode, skip_parameters)
    if engage != None:
        od["parentdir"]["default"] = engage["parentdir"]
        o["parentdir"] = get_run_option("parentdir", od["parentdir"], mode, skip_parameters)
        od["execdir"]["default"] = o["parentdir"]
        od["workdir"]["default"] = o["parentdir"]
        od["histdir"]["default"] = o["parentdir"]
        o["execdir"] = get_run_option("execdir", od["execdir"], mode, skip_parameters)
    else:
        od["parentdir"]["default"] = "."
        o["parentdir"] = get_run_option("parentdir", od["parentdir"], mode, skip_parameters)
        if o["parentdir"] == None:
            temp_mode = "INTERMEDIATE"
            od["execdir"]["default"] = "."
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
    # relative paths typed from here on are relative to the run folder
    root = misc.run_root(o)
    run_dir = os.path.abspath(root) if root else None
    o["modeldir"] = misc.anchor_path(o["modeldir"], run_dir)
    o["utildir"] = os.path.join(o["modeldir"],'scripts')
    od["tgcmdata"]["default"] = tiegcm_env("TIEGCMDATA")
    o["tgcmdata"] = get_run_option("tgcmdata", od["tgcmdata"], mode, skip_parameters, run_dir=run_dir)
    if base_skip == False:
        if benchmark == None and engage == None:
            o["input_file"] = get_run_option("input_file", od["input_file"], mode, skip_parameters,
                                             run_dir=run_dir)
            if o["input_file"]  != None:
                input_build_skip = True
        else:
            o["input_file"] = None
    od["log_file"]["default"] =  f'{o["workdir"]}/{options["simulation"]["job_name"]}.out'
    o["log_file"] = get_run_option("log_file", od["log_file"], mode, skip_parameters, run_dir=run_dir)

    if od["make"]["default"] == None:
        od["make"]["default"] = default_make_file(options["model"]["data"]["modeldir"],
                                                  options["simulation"]["hpc_system"])
    o["make"] = get_run_option("make", od["make"], mode, skip_parameters, run_dir=run_dir)
    od["modelexe"]["default"] = os.path.join(o["execdir"],"tiegcm.exe")
    o["modelexe"] = get_run_option("modelexe", od["modelexe"], mode, skip_parameters, run_dir=run_dir)
    if os.path.isfile(o["modelexe"]) == False:
        o["modelexe"] = in_execdir(o["execdir"], o["modelexe"])
        if args.compile == False and args.onlycompile == False:
            print(f'{YELLOW}Unable to find {o["modelexe"]}, model must be compiled. Use --compile/-c or --onlycompile/-oc {RESET}')
    if args.coupling == True:
        od["coupled_modelexe"]["default"] = os.path.join(o["execdir"],"tiegcm.x")
        o["coupled_modelexe"] = get_run_option("coupled_modelexe", od["coupled_modelexe"], mode, skip_parameters,
                                               run_dir=run_dir)
        if os.path.isfile(o["coupled_modelexe"]) == False:
            o["coupled_modelexe"] = in_execdir(o["execdir"], o["coupled_modelexe"])
            if args.compile == False and args.onlycompile == False:
                print(f'{YELLOW}Unable to find {o["coupled_modelexe"]}, model must be compiled. Use --compile/-c or --onlycompile/-oc {RESET}')

    # every later data-file lookup searches this directory
    if o["tgcmdata"]:
        os.environ["TIEGCMDATA"] = o["tgcmdata"]


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
        # follows the answered mres, not horires' default mres
        o["nres_grid"] = mres_to_nres_grid(o["mres"]) if o.get("mres") is not None else nres_grid

    
    if input_build_skip == True:
        o["segmentation"] = False
    if input_build_skip == False:        

        options["inp"] = {}
        o = options["inp"]

        od = option_descriptions["inp"]     
        od["STEP"]["default"] = STEP

        run_name = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
        # relative to the job's workdir
        histdir = misc.history_dir(options["model"]["data"])
        if benchmark != None:
            oben = benchmarks_options[benchmark]["inp"]
            for on in od:
                if on in oben:
                    if on in ["SOURCE", "GPI_NCFILE", "IMF_NCFILE"] and oben[on] is not None:
                        od[on]["default"] = resolve_benchmark_file(on, benchmark, oben, tiegcm_env("TIEGCMDATA"))
                    elif on in ["OUTPUT", "SECOUT"]:
                        temp_output = oben[on]
                        try:
                            temp_output = temp_output.replace("+histdir+", histdir)
                            temp_output = temp_output.replace("+run_name+", run_name)
                        except:
                            temp_output = temp_output
                        od[on]["default"] = temp_output
                    elif on in ["other_input"]:
                        temp_output = [item.replace("+tiegcmdata+", tiegcm_env("TIEGCMDATA")) if item is not None and item != 'null' else item for item in oben[on]]
                        od[on]["default"] = temp_output
                    elif oben[on] == None:
                        od[on]["default"] = od[on]["default"]
                    else:
                        od[on]["default"] = oben[on]
                
        od["LABEL"]["default"] = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
        if engage is not None:
            # an answer equal to engage's value leaves the key to engage
            for on in COUPLED_DEFAULT_KEYS:
                od[on]["LEVEL"] = "EXPERT"
                od[on]["default"] = engage["coupled_defaults"][on]
        if options["simulation"]["hpc_system"] == "aitken":
            od["SECFLDS"]["warning"] = "Limit SECFLDS. File libraries on aitken are built without big file support."
             
        temp_mode = mode
        skip_inp = []
        segment = None
        start_stop_set = 0
        crosses_year_boundary = False
        mxday_start = 366
        coupled_unset = []                        # coupled default keys answered 'none'
        gswm_defaulted = []
        for on in od:
            if start_stop_set == 0 and benchmark == None:
                temp_mode =  "BASIC"
            if on == "start_time" and benchmark != None:
                continue
            elif on == "stop_time" and benchmark != None:
                continue
            elif on == "start_time" and benchmark == None:
                if engage != None:
                    od["start_time"]["default"] = engage["start_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on == "stop_time" and benchmark == None:
                if engage != None:
                    od["stop_time"]["default"] = engage["stop_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if o["start_time"] is not None and o[on] is not None:
                    o[on] = _reprompt_while(on, od[on], temp_mode, skip_parameters, o[on],
                                            lambda v: _date_order_problem("stop_time", v, o["start_time"], "start_time"))
                if o["start_time"] == None or o[on] == None:
                    raise ValueError(
                        f"{RED}start_time and stop_time are required for a non-benchmark run "
                        f"(got start_time={o['start_time']!r}, stop_time={o[on]!r}).{RESET}"
                    )
                start_stop_set = 1
                temp_mode = mode
                START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(o["start_time"], o["stop_time"])
                od["START_YEAR"]["default"] = START_YEAR
                od["START_DAY"]["default"] = START_DAY
                od["PRISTART"]["default"] = PRISTART
                od["PRISTOP"]["default"] = PRISTOP
                # input.F caps a single run's day-of-year at mxday (366, 367 in a leap year):
                # a run past the year's end must be segmented (split at Jan 1).
                crossing = year_boundary(o["start_time"], o["stop_time"])
                crosses_year_boundary = crossing is not None
                if crossing:
                    mxday_start = crossing[2]
            elif on == "secondary_start_time" and benchmark == None:
                od["secondary_start_time"]["default"] = o["start_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                o[on] = _reprompt_while(on, od[on], temp_mode, skip_parameters, o[on],
                                        lambda v: _date_order_problem("secondary_start_time", v, o["start_time"], "start_time",
                                                                      strictly_after=False, before=o["stop_time"],
                                                                      before_label="stop_time", strictly_before=True))
            elif on == "secondary_stop_time" and benchmark == None:
                od["secondary_stop_time"]["default"] = o["stop_time"]
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                o[on] = _reprompt_while(on, od[on], temp_mode, skip_parameters, o[on],
                                        lambda v: _date_order_problem("secondary_stop_time", v, o["secondary_start_time"],
                                                                      "secondary_start_time", before=o["stop_time"],
                                                                      before_label="stop_time"))
            elif on == "segment" and benchmark == None:
                if engage != None:
                    od["segment"]["default"] = engage["segment"]
                if crosses_year_boundary:
                    yb_warn = (f"This run crosses the {START_YEAR}->{START_YEAR + 1} year boundary "
                               f"(it would run to day {PRISTOP[0]}, but the model caps a single run at "
                               f"day {mxday_start}). Segmentation is required: the run is split "
                               f"automatically at Jan 1. Enter a segment length, e.g. '5 0 0 0'.")
                    od["segment"]["warning"] = (od["segment"]["warning"] + "\n" if od["segment"].get("warning") else "") + yb_warn
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                while crosses_year_boundary and o[on] == [None]:
                    print(f"{RED}Segmentation is required when the run crosses the year boundary "
                          f"(stop day {PRISTOP[0]} > model max {mxday_start} for start year {START_YEAR}). "
                          f"Enter a segment length (Day Hour Min Sec), e.g. '5 0 0 0'.{RESET}")
                    try:
                        o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                    except EOFError:
                        raise ValueError(
                            f"Segmentation is required for a run that crosses the year boundary "
                            f"(stop day {PRISTOP[0]} > model max {mxday_start} for start year "
                            f"{START_YEAR}), but none was provided. Re-run with a segment length "
                            f"(e.g. inp.segment='5 0 0 0')."
                        )
                # segmentation follows the segment length alone
                if parse_segment(o[on]) is None:
                    if options["model"]["specification"].get("segmentation"):
                        print(f"{YELLOW}Segmentation was set to True but no segment length was given: "
                              f"the run is written as a single (unsegmented) job.{RESET}")
                    options["model"]["specification"]["segmentation"] = False
                else:
                    options["model"]["specification"]["segmentation"] = True
                    segment = parse_segment(o[on])
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
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on == "SOURCE":
                if benchmark == None:
                    od["SOURCE"]["default"] = select_source_defaults(options, option_descriptions)    
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if o[on] == None:
                    o[on] = get_run_option(on, od[on], "BASIC", run_dir=run_dir)
            elif on == "SOURCE_START":
                source_mtimes = get_mtime(options["inp"]["SOURCE"])
                if benchmark == None:
                    source_mtimes = _source_start_candidates(o, od, source_mtimes, PRISTART, skip_parameters, run_dir)
                od["SOURCE_START"]["valids"] = source_mtimes
                temp_mode_1 = temp_mode
                if len(od["SOURCE_START"]["valids"]) > 1:    
                    temp_mode_1 = "INTERMEDIATE"
                od["SOURCE_START"]["default"] = od["SOURCE_START"]["valids"][0]
                o[on] = get_run_option(on, od[on], temp_mode_1, skip_parameters, run_dir=run_dir)
                while benchmark == None and _hms_differs(o[on], PRISTART):
                    print(f"{RED}{source_start_hms_error(o[on], PRISTART)}{RESET}")
                    o[on] = get_run_option(on, od[on], temp_mode_1, skip_parameters, run_dir=run_dir)
            elif on == "PRIHIST" and benchmark== None:
                PRIHIST = inp_prihist(PRISTART,PRISTOP, segment, cadence_step(engage, STEP))
                od["PRIHIST"]["default"] =  PRIHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                PRIHIST = [int(i) for i in o[on].split()]
                MXHIST_PRIM_set ,MXHIST_PRIM_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], PRIHIST, od["MXHIST_PRIM"]["warning"], segment)
                od["MXHIST_PRIM"]["default"] = MXHIST_PRIM_set
                od["MXHIST_PRIM"]["warning"] = MXHIST_PRIM_warning_set
            elif on == "MXHIST_PRIM" and benchmark== None:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                MXHIST_PRIM = int(o[on])
            elif on == "OUTPUT" and benchmark== None:
                OUTPUT, pri_files_n = inp_pri_out(o["start_time"], o["stop_time"], PRIHIST, MXHIST_PRIM, 0, histdir,run_name)
                od["OUTPUT"]["default"] = OUTPUT
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on == "SECHIST" and benchmark== None:
                if o.get("secondary_start_time") and o.get("secondary_stop_time"):
                    _, _, SEC_WINDOW_START, SEC_WINDOW_STOP = inp_pri_date(o["secondary_start_time"], o["secondary_stop_time"])
                else:
                    SEC_WINDOW_START, SEC_WINDOW_STOP = PRISTART, PRISTOP
                SECHIST = inp_sechist(SEC_WINDOW_START, SEC_WINDOW_STOP, segment, cadence_step(engage, o.get("STEP"), STEP))
                od["SECHIST"]["default"] =  SECHIST
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                SECHIST = [int(i) for i in o[on].split()]
                MXHIST_SECH_set ,MXHIST_SECH_warning_set  = inp_mxhist(o["start_time"], o["stop_time"], SECHIST, od["MXHIST_SECH"]["warning"],segment)
                od["MXHIST_SECH"]["default"] = MXHIST_SECH_set
                od["MXHIST_SECH"]["warning"] = MXHIST_SECH_warning_set
                SECSTART, SECSTOP = inp_sec_date(o["secondary_start_time"], o["secondary_stop_time"], SECHIST)
                od["SECSTART"]["default"] = SECSTART
                od["SECSTOP"]["default"] = SECSTOP                
            elif on == "MXHIST_SECH" and benchmark== None:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                MXHIST_SECH = int(o[on])
            elif on == "SECOUT" and benchmark== None:
                SECOUT, sec_files_n = inp_sec_out(o["secondary_start_time"], o["secondary_stop_time"],  SECHIST, MXHIST_SECH, 0, histdir,run_name)
                od["SECOUT"]["default"] = SECOUT
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on == "POTENTIAL_MODEL":
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if o[on] == "HEELIS":
                    skip_inp_temp = ["IMF_NCFILE","BXIMF","BYIMF","BZIMF","SWDEN","SWVEL"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
                elif o[on] == "WEIMER":
                    # Weimer computes CTPOTEN, and input.F refuses KP under WEIMER
                    skip_inp_temp = ["CTPOTEN", "KP"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on == "ONEWAY":
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)            
            elif on == "GPI_NCFILE" and on not in skip_inp:
                if benchmark == None and od["GPI_NCFILE"]["default"] is None:
                    # newest bundled GPI if it covers the run window, else generate one
                    latest = select_latest_gpi(tiegcm_env("TIEGCMDATA"))
                    if latest is not None and file_covers(latest, o["start_time"], o["stop_time"]):
                        od["GPI_NCFILE"]["default"] = latest
                    else:
                        od["GPI_NCFILE"]["default"] = "gen"
                    gen_note = ("Enter 'gen' to generate a GPI file for this run's dates "
                                "with gcmprocpy (27-day trailing F10.7 average).")
                    od["GPI_NCFILE"]["warning"] = (od["GPI_NCFILE"]["warning"] + "\n" if od["GPI_NCFILE"].get("warning") else "") + gen_note
                gpi_generated = False
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                while isinstance(o[on], str) and o[on].strip().lower() == "gen":
                    try:
                        o[on] = generate_gpi_file(o["start_time"], o["stop_time"], options["model"]["data"]["workdir"])
                        gpi_generated = True
                    except Exception as e:
                        print(f"{RED}GPI generation failed: {e}{RESET}")
                        print(f"{YELLOW}Enter a GPI file path instead, or 'none' to skip GPI.{RESET}")
                        od["GPI_NCFILE"]["default"] = None  # don't re-trigger gen on ENTER
                        o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if gpi_generated and isinstance(o[on], str) and not file_covers(o[on], o["start_time"], o["stop_time"]):
                    print(f"{YELLOW}Warning: generated GPI {o[on]} does not span the full run "
                          f"({o['start_time']}..{o['stop_time']}); the indices source may lag real time.{RESET}")
                if o[on] != None:
                    skip_inp_temp = ["KP","POWER","CTPOTEN","F107","F107A"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
                    od["F107"]["warning"] = "F10.7 can be read by GPI File and can be skipped."
                    od["F107A"]["warning"] = "81-Day Average of F10.7 can be read by GPI File and can be skipped."
            elif on == "IMF_NCFILE" and on not in skip_inp:
                if benchmark == None:
                    gen_note = ("Enter 'gen' to generate an IMF file for this run's dates "
                                "with gcmprocpy (OMNI solar wind).")
                    od["IMF_NCFILE"]["warning"] = (od["IMF_NCFILE"]["warning"] + "\n" if od["IMF_NCFILE"].get("warning") else "") + gen_note
                imf_generated = False
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                while isinstance(o[on], str) and o[on].strip().lower() == "gen":
                    try:
                        o[on] = generate_imf_file(o["start_time"], o["stop_time"], options["model"]["data"]["workdir"])
                        imf_generated = True
                    except Exception as e:
                        print(f"{RED}IMF generation failed: {e}{RESET}")
                        print(f"{YELLOW}Enter an IMF file path instead, or 'none' to skip IMF.{RESET}")
                        od["IMF_NCFILE"]["default"] = None  # don't re-trigger gen on ENTER
                        o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if imf_generated and isinstance(o[on], str) and not file_covers(o[on], o["start_time"], o["stop_time"]):
                    print(f"{YELLOW}Warning: generated IMF {o[on]} does not span the full run "
                          f"({o['start_time']}..{o['stop_time']}); the OMNI source may lag real time.{RESET}")
                if o[on] != None:
                    skip_inp_temp = ["BXIMF","BYIMF","BZIMF","SWDEN","SWVEL"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on == "KP" and on not in skip_inp:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
                if o[on] != None:
                    skip_inp_temp = ["POWER","CTPOTEN"]
                    for item in skip_inp_temp:
                        if item not in skip_inp:
                            skip_inp.append(item)
            elif on in GSWM_PATTERNS:
                if not gswm_allowed(on, horires):
                    o[on] = None
                    continue
                gswm_file = find_file(GSWM_PATTERNS[on].format(horires=horires), tiegcm_env("TIEGCMDATA"))
                od[on]["default"] = f"{gswm_file}" if gswm_file is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on == "HE_COEFS_NCFILE":
                # a missing file stays None, never the string 'None'
                HE_COEFS_NCFILE = find_file(he_coefs_pattern(horires), tiegcm_env("TIEGCMDATA"))
                od[on]["default"] = f"{HE_COEFS_NCFILE}" if HE_COEFS_NCFILE is not None else None
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on not in skip_inp:
                o[on] = get_run_option(on, od[on], temp_mode, skip_parameters, run_dir=run_dir)
            elif on in skip_inp:
                o[on] = od[on]["default"]
            if engage is not None and on in COUPLED_DEFAULT_KEYS:
                if str(o[on]) == str(od[on]["default"]):
                    o[on] = None                  # left to engage
                elif o[on] is None:
                    coupled_unset.append(on)      # answered 'none': stays unset
            elif on in GSWM_PATTERNS and o[on] is not None and str(o[on]) == str(od[on]["default"]):
                gswm_defaulted.append(on)
            if on in DERIVED_INP_KEYS and _prompted(on, od[on], temp_mode, skip_parameters):
                segmented = engage is not None or bool(options["model"]["specification"].get("segmentation"))
                if _derived_note(on, o[on], od[on]["default"], segmented):
                    o[on] = get_run_option(on, od[on], "BASIC")      # the derived value
        _complete_forcing(o, od, "BASIC" if temp_mode == "BENCH" else temp_mode, skip_parameters)
        if engage is not None:
            engage["coupled_unset"] = coupled_unset
            apply_coupled_defaults(o, horires, coupled_unset)
        elif lbc_other_set(o) and any(o.get(k) is not None for k in gswm_defaulted):
            # default GSWM files give way to another lower-boundary source
            for k in gswm_defaulted:
                o[k] = None
            print(f"{YELLOW}NOTE: the default GSWM files are left out: another lower-boundary "
                  f"source is set.{RESET}")


    hpc_platform = options["simulation"]["hpc_system"]
    if hpc_platform == "linux":
        pbs_build_skip = True
    if pbs_build_skip == False:
        options["job"] = {}
        o = options["job"]
        skip_pbs = []
        hpc_platform = options["simulation"]["hpc_system"]
        od = option_descriptions["job"][hpc_platform]
        if engage != None:
            for key in ("queue", "job_priority", "walltime", "group_list"):
                if key in od and key in engage:
                    od[key]['default'] = engage[key]
        for on in od:
            if on == "resource":
                options["job"]["resource"] = {}
                odt = od["resource"]
                ot = options["job"]["resource"]
                if engage != None and engage.get("model") and "model" in odt:
                    odt["model"]["default"] = engage["model"]
                for ont in odt:
                    if hpc_platform != "aitken":
                        select_default,ncpus_default,mpiprocs_default = select_resource_defaults(options,option_descriptions)
                        odt["select"]["default"] = select_default
                        odt["ncpus"]["default"] = ncpus_default
                        odt["mpiprocs"]["default"] = mpiprocs_default
                    else:
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
                        # no more ranks per node than the cores asked for
                        if ot.get("ncpus") is not None and odt[ont]["default"] not in (None, "ncpus"):
                            odt[ont]["default"] = min(int(odt[ont]["default"]), int(ot["ncpus"]))
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
                        mpiprocs = ot[ont]
                    elif ont != "model" or "model" not in ot:
                        # aitken's model was asked above; the select/ncpus defaults derive from it
                        ot[ont] = get_run_option(ont, odt[ont], mode, skip_parameters)
            elif on =="nprocs":
                od[on]["default"] = int(nnodes) * int(mpiprocs)
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
            elif on == "project_code":
                if engage != None:
                    od[on]["default"] = engage["project_code"]
                    od[on].pop("choices", None)
                elif od[on].get("default") in (None, "") and len(od[on].get("choices") or []) == 1:
                    od[on]["default"] = od[on]["choices"][0]       # the user's only project code
                o[on] = _reprompt_while(
                    on, od[on], mode, skip_parameters, get_run_option(on, od[on], mode, skip_parameters),
                    lambda v: None if str(v).strip() not in ("", "None", "none", "null") else
                    (f"job.project_code is required on {hpc_platform}: enter the PBS account (#PBS -A) "
                     f"the jobs are charged to" + (" (engage's pbs.account_name is not set)" if engage else "")))
            elif on not in skip_pbs:
                o[on] = get_run_option(on, od[on], mode, skip_parameters)
                if on == "queue" and engage is None and "walltime" in od and o[on]:
                    od["walltime"]["default"] = (queue_walltime_default(hpc_platform, o[on])
                                                 or od["walltime"]["default"])

    return options

def regrid_source(options, run_name):
    """Regrid inp.SOURCE onto the run's grid as <workdir>/<run_name>_prim.nc; return that path."""
    out_prim = source_target(options, run_name)
    regrid(options, out_prim)
    return out_prim


def regrid(options, out_prim):
    """Regrid inp.SOURCE onto the run's grid (horires, vertres, zitop) as out_prim."""
    spec = options["model"]["specification"]
    interpic(options["inp"]["SOURCE"], float(spec["horires"]), float(spec["vertres"]),
             float(spec["zitop"]), out_prim)


def source_target(options, run_name):
    """Return the regridded SOURCE path <workdir>/<run_name>_prim.nc (written by regrid).

    Always regridded: an existing file may come from an earlier run with another SOURCE or zitop.
    """
    workdir = options["model"]["data"]["workdir"]
    source = options["inp"]["SOURCE"]
    out_prim = f'{workdir}/{run_name}_prim.nc'
    # SOURCE_START must be an exact history of SOURCE, else the model stops.
    try:
        source_mtimes = get_mtime(source)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError):
        source_mtimes = None        # an unreadable SOURCE fails in interpic with its own error
    if source_mtimes is not None:
        value, note = resolve_source_start(source, options["inp"].get("SOURCE_START"), source_mtimes)
        if note:
            print(f"{YELLOW}{note}{RESET}")
        options["inp"]["SOURCE_START"] = value
    if source is not None and os.path.abspath(source) == os.path.abspath(out_prim):
        raise ValueError(f"{RED}inp.SOURCE {source} is this run's own regridded initial condition, "
                         f"which the regrid would overwrite: set inp.SOURCE to the original "
                         f"initial-condition file.{RESET}")
    return out_prim

def _exe_path(options, key, invoke_dir):
    exe = options.get("model", {}).get("data", {}).get(key)
    if exe in (None, "", "None"):
        return None
    return os.path.normpath(os.path.join(invoke_dir, str(exe)))


def _exe_grids(options, args):
    """Return {role: (coupled, horires, vertres, nres_grid)} for the executables the run launches.

    Under engage the coupled executable is at engage's horires_coupled grid.
    """
    spec = options["model"]["specification"]
    grids = {"modelexe": (False, spec.get("horires"), spec.get("vertres"), spec.get("nres_grid"))}
    if args.coupling:
        if isinstance(args.engage, dict) and args.engage.get("horires_coupled") is not None:
            h = float(args.engage["horires_coupled"])
            v, _, n, _ = resolution_solver(h)
            grids["coupled_modelexe"] = (True, h, v, n)
        else:
            grids["coupled_modelexe"] = (True, spec.get("horires"), spec.get("vertres"), spec.get("nres_grid"))
    return grids


def check_executables(options, args, invoke_dir):
    """Check the run's executables against the run's build kind and grid; return the warnings.

    A non-ELF file is refused; a wrong kind or grid is refused under engage, a warning otherwise.
    """
    if args.onlycompile:
        return []
    building = ("coupled_modelexe" if args.coupling else "modelexe") if args.compile else None
    zitop = options["model"]["specification"].get("zitop")
    engage = args.engage is not None
    refused, warned = [], []
    for role, (coupled, horires, vertres, nres_grid) in _exe_grids(options, args).items():
        path = _exe_path(options, role, invoke_dir)
        if role == building or path is None:
            continue
        not_elf, problems = misc.exe_problems(misc.exe_identity(path, zitop), f"model.data.{role}",
                                              coupled, horires, vertres, nres_grid, zitop)
        if not_elf:
            refused.append(not_elf)
        if engage:
            refused.extend(problems)
        elif role == "coupled_modelexe":
            warned.extend(f"{p} (this standalone job does not launch it; an engage run would refuse it)"
                          for p in problems)
        else:
            warned.extend(problems)
    if refused:
        raise ValueError(f"{RED}The TIE-GCM executable cannot run this "
                         f"{'coupled ' if engage else ''}run:\n  - " + "\n  - ".join(refused)
                         + f"\nPoint the option at the right executable, or build it (--compile / "
                           f"--onlycompile{' --coupling' if args.coupling else ''}).{RESET}")
    for w in warned:
        print(f"{YELLOW}WARNING: {w}.{RESET}")
    return warned


def tiegcmrun(args=None):
    parser = create_command_line_parser()
    if args is not None:
        args = parser.parse_args(args)
    else:
        args = parser.parse_args()
    args.tui = use_tui(args.no_tui)
    rederive = ()
    if args.rederive is not None:
        try:
            rederive = parse_rederive(args.rederive)
        except ValueError as e:
            parser.error(str(e))
        if rederive == ("list",):
            print(rederive_groups_text())
            return None
        if args.options_path is None and args.engage is None:
            parser.error("--rederive needs -o <options JSON>")
    _require_tiegcm_env()
    # relative paths in the options resolve against this, even after the chdir into workdir
    invoke_dir = os.getcwd()
    if args.engage is not None:
        args.coupling = True
    clobber = args.clobber
    debug = args.debug
    debug_build = args.debug_build
    options_path = args.options_path
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
    if benchmark != None and mode == None:
        mode = "BENCH"
    elif mode == None:
        mode = "BASIC"
    compile_flag = compile or onlycompile
    submit_all_jobs_script = ''
    pbs_script = None
    linux_run_script = None
    linux_hints = None
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
        # relative paths in the JSON are relative to its own directory
        options = read_options_json(options_path)
        problem = options_file_problem(options)
        if problem:
            raise ValueError(f"{RED}{options_path}: {problem}{RESET}")
        misc.resolve_json_paths(options, os.path.dirname(os.path.abspath(options_path)))
        coupled_unset = [k for k in (options.pop("coupled_unset", None) or ()) if k in COUPLED_DEFAULT_KEYS]
        # a typed "none"/"null" string is an explicit 'none'; a JSON null is not
        inp = options.get("inp") if engage is not None and isinstance(options.get("inp"), dict) else {}
        for k in COUPLED_DEFAULT_KEYS:
            if isinstance(inp.get(k), str) and misc.is_unset(inp[k]):
                inp[k] = None
                if k not in coupled_unset:
                    coupled_unset.append(k)
        with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
            option_descriptions = apply_system_config(json.load(f), socket.gethostname())
        if engage != None:
            # engage owns the run window
            rederive_pop(options, [g for g in rederive if g != "dates"])
            args.engage["coupled_unset"] = coupled_unset
            options = engage_options_updater(options, args.engage, option_descriptions)
        else:
            options = prepare_replay(options, option_descriptions, rederive)
    else:
        options = None
        if getattr(args, "tui", False):
            _tf = _import_own_tui_form()
            try:
                options = _tf.run_tui_wizard(args)
            except _tf._FormCancelled:
                print(f"{YELLOW}tiegcmrun wizard cancelled (Ctrl-C); tiegcmrun wrote nothing.{RESET}")
                raise SystemExit(130) from None
            except _tf._FormUnavailable as e:
                print(f"{YELLOW}Wizard unavailable ({e}); using the linear prompts.{RESET}")
                if debug:
                    import traceback
                    traceback.print_exc()
            if options is not None:
                # side effects the linear prompts do inline (GPI/IMF "gen", TIEGCMDATA, exe check)
                options = _tf.apply_accept_actions(options, args)
        if options is None:
            options = prompt_user_for_run_options(args)
    generated = COUPLED_HISTORY_KEYS if engage is not None else ()
    if args.check:
        inp = options.get("inp") if isinstance(options.get("inp"), dict) else {}
        pending = tuple(f"inp.{k}" for k in ("GPI_NCFILE", "IMF_NCFILE")
                        if isinstance(inp.get(k), str) and inp[k].strip().lower() == "gen")
        validate(options, generated=generated + pending)
        check_executables(options, args, invoke_dir)
        return (options, [], [])
    # a replayed GPI/IMF 'gen' is generated here
    if options_path and isinstance(options.get("inp"), dict):
        misc.generate_pending(options["inp"], options["model"]["data"].get("workdir"))
    validate(options, generated=generated)
    if debug:
        print(f"options = {options}")

    check_executables(options, args, invoke_dir)

    run_name = f"{options['simulation']['job_name']}_{options['model']['specification']['horires']}x{options['model']['specification']['vertres']}"
    execdir   = options["model"]["data"]["execdir"]
    workdir = options["model"]["data"]["workdir"]
    outdir = options["model"]["data"]["histdir"]
    for d in (workdir, outdir, execdir):
        os.makedirs(d, exist_ok=True)
    workdir_abs = os.path.abspath(workdir)
    json_path = f"{workdir}/{PARAMETERS_JSON}"
    # checked before anything is written; engage gates its own files
    if args.onlycompile == False and args.engage == None and os.path.exists(json_path) and not clobber:
        raise FileExistsError(f"Options file {json_path} exists! Use --clobber to overwrite it and "
                              f"the run's .inp/.pbs files.")

    if options.get("inp") == None:
        input_file_generatred = True
    else:
        input_file_generatred = False
    if args.onlycompile == True:
        compile_tiegcm(options, debug_build, coupling, hidra)
    elif args.engage != None:
        options_coupling,standalone_pbs_files,coupling_inp_files = engage_run(options, debug, coupling, args.engage)
        return (options_coupling,standalone_pbs_files,coupling_inp_files)
    else:
        if args.compile == True:
            compile_tiegcm(options, debug_build, coupling, hidra)
        # segmentation follows inp.segment alone
        inp_options = options.get("inp")
        options["model"]["specification"]["segmentation"] = (
            bool(inp_options) and parse_segment(inp_options.get("segment")) is not None)
        if options["model"]["specification"]["segmentation"] == False:
            # The job is rendered before any file is written, so a refused job leaves nothing.
            if input_file_generatred == False:
                regridded = source_target(options, run_name)
                options["model"]["data"]["input_file"] = run_file(workdir, run_name, None, "inp")
            if options["model"]["data"]["log_file"] == None:
                options["model"]["data"]["log_file"] = os.path.join( options["model"]["data"]["workdir"], f"{run_name}.out")
            pbs_text = render_pbs(options) if options["simulation"]["hpc_system"] != "linux" else None
            if input_file_generatred == False:
                # the saved options keep the user's inp.SOURCE
                regrid(options, regridded)
                inp_options = copy.deepcopy(options)
                inp_options["inp"]["SOURCE"] = regridded
                create_inp_scripts(inp_options,run_name,None)

            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(misc.json_paths(options, workdir_abs), f, indent=JSON_INDENT)
            print(saved_parameters_note(workdir_abs, invoke_dir))

            if pbs_text is not None:
                pbs_script = create_pbs_scripts(options, run_name, None, content=pbs_text)
        else:
            first_source = source_target(options, run_name)
            planned = segment_jobs(options, run_name, pbs=True, first_source=first_source)
            regrid(options, first_source)
            with open(json_path, "w", encoding="utf-8") as f:
                json.dump(misc.json_paths(options, workdir_abs), f, indent=JSON_INDENT)
            print(saved_parameters_note(workdir_abs, invoke_dir))
            inp_files, pbs_files,log_files, pristart_times, pristop_times, _ = write_segment_jobs(
                planned, run_name)
            init_inp = inp_files[0]    
            init_pbs = pbs_files[0]
            options["model"]["data"]["input_file"] = init_inp

            if options["simulation"]["hpc_system"] == "linux":
                linux_commands = linux_run_commands(options, inp_files, log_files, relative=True)
                linux_hints = linux_run_commands(options, inp_files, log_files)
                linux_run_script = os.path.join(workdir, f"{options['simulation']['job_name']}_run.sh")
                with open(linux_run_script, "w", encoding="utf-8") as f:
                    f.write("\n".join(["#!/bin/bash", "set -euo pipefail", cdpath_reset("bash"),
                                       'cd "$(dirname "$0")"'] + linux_commands) + "\n")
                os.chmod(linux_run_script, 0o755)
            else:
                # stops at the first rejected submission
                os.chdir(workdir_abs)
                submit_all_jobs_script = (f"{options['simulation']['job_name']}_pbs.sh")
                jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"],
                            node_type=options["simulation"].get("node_type"))
                chain = [(".", os.path.basename(p)) for p in pbs_files if p]
                with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
                    f.write(jg.submit_chain_script(
                        [("segments", "tiegcm_job_id", chain, [])],
                        header=f"Submit the {options['simulation']['job_name']} TIE-GCM run. Run from anywhere.",
                        jobids=f"{options['simulation']['job_name']}.jobids"))
                os.chmod(submit_all_jobs_script, 0o755)

    if args.execute == True and args.onlycompile == False and options["simulation"]["hpc_system"] != "linux":
        jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"],
                    node_type=options["simulation"].get("node_type"))
        # a --coupling build produces tiegcm.x, not the standalone modelexe
        if not os.path.exists(os.path.join(invoke_dir, options["model"]["data"]["modelexe"])):
            print(f'{RED}Unable to find executable {options["model"]["data"]["modelexe"]} in {execdir}{RESET}')
            if args.coupling:
                print(f'{YELLOW}--coupling builds the coupled executable, which a standalone job '
                      f'cannot run; build without --coupling to submit a standalone run.{RESET}')
            exit(1)
        # submitted from the workdir so $PBS_O_WORKDIR is the job's directory
        submit_env = dict(os.environ, PWD=workdir_abs)
        try:
            if submit_all_jobs_script == '':
                result = subprocess.run(jg.b["submit"].split() + ["./" + os.path.basename(pbs_script)], check=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                        cwd=workdir_abs, env=submit_env)
                job_id = result.stdout.strip()
                print(f'Job submitted successfully. Job ID: {job_id}')
            else:
                result = subprocess.run(['./'+submit_all_jobs_script], check=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True, cwd=workdir_abs, env=submit_env)
                job_ids = " ".join(result.stdout.split())
                print(f'Jobs submitted successfully. Job IDs: {job_ids}')
        except subprocess.CalledProcessError as e:
            print(f'{YELLOW}Error submitting job: {e.stderr}{RESET}')
            print(f"{YELLOW}Check the job script for errors{RESET}")
            print(f"To submit job use command {YELLOW}{_submit_hint(jg, pbs_script, submit_all_jobs_script, workdir_abs, invoke_dir)}{RESET}")

    elif args.onlycompile == False and options["simulation"]["hpc_system"] != "linux":
        jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"],
                    node_type=options["simulation"].get("node_type"))
        print(f"{YELLOW}Execute is set to false{RESET}")
        print(f"To submit job use command {YELLOW}{_submit_hint(jg, pbs_script, submit_all_jobs_script, workdir_abs, invoke_dir)}{RESET}")
    elif args.onlycompile == False:
        print(f"{YELLOW}HPC System is set to linux{RESET}")
        if linux_run_script:
            print(f"{YELLOW}To run every segment in order use command{RESET} {linux_run_script}")
            print(f"{YELLOW}or run the segments one at a time, in this order (from {workdir_abs}):{RESET}")
            linux_commands = linux_hints
        else:
            print(f"{YELLOW}From {workdir_abs}:{RESET}")
            linux_commands = linux_run_commands(options, [options['model']['data']['input_file']],
                                                [options['model']['data']['log_file']])
        for command in linux_commands:
            print(f"{YELLOW}To run the model use command{RESET} {command}")


def options_file_problem(options, unknown=False):
    """Return why a loaded options JSON is refused, or None; `unknown` also refuses unknown keys."""
    problem = section_problem(options) or old_format_problem(options)
    if problem or not unknown:
        return problem
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = apply_system_config(json.load(f), socket.gethostname())
    return unknown_keys_problem(options, options_key_schema(option_descriptions))


def validate(options, generated=()):
    """Raise ValueError listing every problem in the options (`generated`: keys the run derives)."""
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = apply_system_config(json.load(f), socket.gethostname())
    problems = validate_options(options, option_descriptions, generated=generated)
    if len(problems) == 1:
        raise ValueError(f"{RED}{problems[0]}{RESET}")
    if problems:
        raise ValueError(f"{RED}The options cannot be run:\n  - " + "\n  - ".join(problems) + RESET)
    return options


def saved_parameters_note(workdir_abs, invoke_dir):
    """Return the line naming the saved options JSON and its replay command."""
    saved = spell(os.path.join(workdir_abs, PARAMETERS_JSON), invoke_dir, invoke_dir)
    return f"The options are saved in {saved}; replay them with tiegcmrun.py -o {saved}."


def _submit_hint(jg, pbs_script, submit_all_jobs_script, workdir_abs, invoke_dir):
    """Return the submit command, runnable from the directory tiegcmrun ran in."""
    if submit_all_jobs_script == '':
        return jg.submit_command(os.path.basename(pbs_script),
                                 subdir=spell(workdir_abs, invoke_dir, invoke_dir))
    return spell(os.path.join(workdir_abs, submit_all_jobs_script), invoke_dir, invoke_dir)

def main(argv=None):
    """Run tiegcmrun(argv); print a refusal as one ERROR line (traceback only with --debug)."""
    argv = sys.argv[1:] if argv is None else list(argv)
    try:
        tiegcmrun(argv)
    except (ValueError, FileExistsError, RuntimeError) as e:
        if "--debug" in argv or "-d" in argv:
            raise
        text = re.sub(r"\x1b\[[0-9;]*m", "", str(e))
        print(f"{RED}ERROR: {text}{RESET}")
        return 1
    return 0


if __name__ == "__main__":
    """Begin tiegcmrun program."""
    sys.exit(main())
