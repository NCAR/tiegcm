"""
TIE-GCM side of a coupled engage (GRT) run: reads engage's options and writes the TIE-GCM
spin-up jobs and the coupled namelists.
"""

import os
import json
import copy
import math
import socket
from datetime import datetime, timedelta

from misc import (seconds_to_dhms, resolution_solver, mres_to_nres_grid, find_file, select_resource_defaults,
                  get_mtime, select_source_defaults, segment_time, gswm_allowed, ntask_error, json_paths,
                  tiegcm_env, is_unset, GSWM_PATTERNS)
from output_solver import (segment_jobs, write_segment_jobs, segment_plan, coupled_handoff_regrids, CONFIG_DIR,
                           job_resource_problems, apply_system_config)
from interpolation import interpic
from namelist_solver import (cadence_problems, inp_pri_date, nc_wrhist_mtime, parse_segment, parse_run_datetime, resolve_source_start,
                             lbc_other_set)
from replay_solver import (fill_engage_omitted, secflds_replay, machine_hint, names_hint, validate_options,
                           ensure_sections)
from jobgen import JobGen

YELLOW = '\033[33m'
RESET = '\033[0m'

RESOURCE_KEYS = ("select", "ncpus", "mpiprocs")

# Coupled-run defaults for the O+ / temperature caps (prompted only at EXPERT).
COUPLED_CAPS = {"OPDIFFCAP": "2e9", "OPDIFFRATE": "0.3", "OPDIFFLEV": "7", "OPFLOOR": "3000",
                "OPRATE": "0.3", "OPLEV": "7", "OPLATWIDTH": "20", "TE_CAP": "8000", "TI_CAP": "8000"}
COUPLED_DEFAULT_KEYS = tuple(COUPLED_CAPS) + tuple(GSWM_PATTERNS)

OPTION_DESCRIPTIONS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "options_description.json")

def gamres_to_res(gamres):
    # GAMERA grid -> (spin-up horires, coupled horires) in degrees; H shares O's since 0.625 is
    # TIE-GCM's finest practical grid.
    if gamres == "D":
        return 2.5 , 2.5
    elif gamres == "Q":
        return 2.5 , 1.25
    elif gamres == "O":
        return 1.25 , 0.625
    elif gamres == "H":
        return 1.25 ,0.625
    else:
        raise ValueError(f"Unknown gamera_grid_type: {gamres!r} (expected D, Q, O, or H)")

def engage_parser(engage_parameters):
    """Flatten engage's nested options into the dict tiegcmrun uses; checks the coupled cadence."""

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
    account_name = o.get('account_name')
    # machines without a scheduler have no queue / walltime
    queue = o.get('queue')
    if hpc_system == "derecho":
        job_priority = o['job_priority']
    elif hpc_system == 'aitken':
        group_list = o['group_list']
        node_type = engage_node_model(hpc_system, engage_parameters["simulation"].get("node_type"))
    walltime = o.get('walltime')
    run_directory = o.get('run_directory')

    o = engage_parameters["coupling"]

    # seconds may be stored as '14400.0'
    gr_warm_up_time = int(float(o['gr_warm_up_time']))
    gcm_spin_up_time = int(float(o['gcm_spin_up_time']))
    conda_env = o['conda_env']
    
    start_date = get_engage_start_time(coupled_start_date,gr_warm_up_time+gcm_spin_up_time)
    
    o = engage_parameters["voltron"]
    STEP, voltron_dtOut = coupled_cadence(o['coupling']['dtCouple'], o['output']['dtOut'],
                                          segment_duration, coupled_start_date, stop_date)

    root_directory = os.path.abspath(run_directory if run_directory else os.curdir)
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

    eo['project_code'] = account_name
    if hpc_system == 'aitken':
        eo['group_list'] = group_list
    else:
        eo['group_list'] = None
    eo['queue'] = queue
    
    eo['walltime'] = walltime
    eo['conda_env'] = conda_env
    
    if hpc_system == "derecho":
        eo['job_priority'] = job_priority
    elif hpc_system == 'aitken':
        eo['model'] = node_type
    
    eo["coupled_defaults"] = coupled_inp_defaults({}, horires_standalone)
    eo['skip']= ['group_list','job_name','hpc_system','horires','parentdir','vertres', 'mres', 'input_file', 'LABEL','start_time','stop_time','secondary_start_time','secondary_stop_time','segment' ,'SOURCE_START','PRIHIST','MXHIST_PRIM','SECHIST','MXHIST_SECH','project_code','queue','job_priority','model','walltime','ONEWAY']
    
    return engage_options

def _whole_seconds(value, name):
    """A voltron time ('5.0', 60, ...) as int seconds; ValueError if not whole."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"voltron {name} must be a number of seconds; got {value!r}") from None
    if not f.is_integer():
        raise ValueError(f"voltron {name} = {value} s is not a whole number of seconds: TIE-GCM's "
                         f"STEP is a whole number, so a fractional value would make the TIE-GCM and "
                         f"voltron clocks drift apart. Use whole seconds.")
    return int(f)


def coupled_cadence(dtCouple, dtOut, segment_seconds, coupled_start, stop):
    """Return (STEP, SECHIST seconds) = (dtCouple, dtOut) after checking them against the
    TIE-GCM history rules in input.F; ValueError otherwise."""
    step = _whole_seconds(dtCouple, "coupling.dtCouple")
    if step <= 0 or 60 % step != 0:
        raise ValueError(
            f"voltron coupling.dtCouple = {dtCouple} s cannot be the coupled TIE-GCM STEP: it must "
            f"be a whole number of seconds that divides 60 (1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30 or "
            f"60), because every TIE-GCM history time and cadence must be a multiple of STEP.")
    out = _whole_seconds(dtOut, "output.dtOut")
    if out <= 0 or out % step != 0:
        raise ValueError(
            f"voltron output.dtOut = {dtOut} s is not a multiple of the coupled STEP (dtCouple) "
            f"{step} s: it is the coupled TIE-GCM secondary-history cadence (SECHIST), which must be "
            f"a multiple of STEP. Use a multiple of {step} s.")
    segment_seconds = int(segment_seconds)
    if segment_seconds <= 0 or segment_seconds % out != 0:
        raise ValueError(
            f"The coupled segment_duration {segment_seconds} s is not a multiple of voltron "
            f"output.dtOut {out} s: each coupled segment writes secondary histories every dtOut from "
            f"its start to its stop, so dtOut must divide the segment.")
    for name, value in (("start_date", coupled_start), ("stop_date", stop)):
        t = parse_run_datetime(value, f"simulation.{name}")
        second_of_day = t.hour * 3600 + t.minute * 60 + t.second
        if second_of_day % step != 0:
            raise ValueError(
                f"simulation.{name} = {value} is not on the coupled STEP clock: its time of day "
                f"({second_of_day} s) is not a multiple of dtCouple {step} s. Move it to a "
                f"multiple of {step} s.")
    return step, out


def engage_node_model(hpc_system, node_type=None):
    """The PBS node :model for `node_type` (or the machine default); None without node types."""
    try:
        jg = JobGen(CONFIG_DIR, machine=hpc_system, node_type=node_type or None)
    except ValueError as e:
        raise ValueError(f"engage simulation.node_type={node_type!r} on {hpc_system}: {e} "
                         f"(tiegcmrun config/machines.yaml)") from None
    return jg._node_model()


# Set per segment by engage_run, so not validated from a -to JSON.
COUPLED_HISTORY_KEYS = ("inp.PRIHIST", "inp.MXHIST_PRIM", "inp.SECHIST", "inp.MXHIST_SECH")


def engage_options_updater(options, engage_options, option_descriptions):
    """Apply engage-owned values to a replayed TIE-GCM options JSON (engage -to) and fill omitted
    keys; returns the options."""
    ensure_sections(options)
    saved = copy.deepcopy(options)
    try:
        secflds_replay(options, where="inp.SECFLDS (TIE-GCM options JSON)")
    except ValueError as e:
        raise ValueError(f"\033[31m{e}{RESET}") from None
    o = options["simulation"]
    o["job_name"] = engage_options["job_name"]
    o["hpc_system"] = engage_options["hpc_system"]
    hint = machine_hint(options, option_descriptions)
    if hint:
        print(f"{YELLOW}{hint}{RESET}")
    o = options["model"]["data"]
    o["parentdir"] = engage_options["parentdir"]
    o["execdir"] = o["parentdir"]
    o["workdir"] = o["parentdir"]
    o["histdir"] = o["parentdir"]
    if isinstance(o.get("log_file"), str) and o["log_file"].strip() not in ("", "None", "none", "null"):
        o["log_file"] = os.path.join(o["parentdir"], os.path.basename(o["log_file"]))
    o = options["model"]["specification"]
    o["horires"] = engage_options["horires"]
    vertres, mres, nres_grid, STEP = resolution_solver(o["horires"])
    if o.get("vertres") is None:
        o["vertres"] = vertres
    if o.get("mres") is None:
        o["mres"] = mres
    if o.get("nres_grid") is None:
        o["nres_grid"] = mres_to_nres_grid(o["mres"])
    fill_engage_omitted(options, option_descriptions, engage_options["skip"], engage_options)
    apply_coupled_defaults(options["inp"], engage_options["horires"],
                           engage_options.get("coupled_unset", ()))
    hint = names_hint(options, label_only=True)
    if hint:
        print(f"{YELLOW}{hint}{RESET}")
    o = options["inp"]
    if o.get("STEP") is None:
        o["STEP"] = STEP
    o["start_time"] = engage_options["start_time"]
    o["stop_time"] = engage_options["stop_time"]
    o["secondary_start_time"] = engage_options["start_time"]
    o["secondary_stop_time"] = engage_options["stop_time"]
    o["segment"] = " ".join(map(str, engage_options["segment"]))
    # history cadences are derived per segment in engage_run; only note a given value that won't fit
    try:
        unfit = cadence_problems({k: o.get(k) for k in ("PRIHIST", "SECHIST")}, parse_segment(o["segment"]),
                                 o["STEP"])
    except (ValueError, TypeError):
        unfit = ["inp.PRIHIST", "inp.SECHIST"]
    given = [f"inp.{k} {o[k]}" for k in ("PRIHIST", "SECHIST")
             if o.get(k) is not None and any(p.startswith(f"inp.{k}") for p in unfit)]
    if given:
        print(f"{YELLOW}NOTE: {', '.join(given)} not used: a coupled run derives "
              f"{'it' if len(given) == 1 else 'them'} per segment.{RESET}")
    # one-way coupling would deadlock the voltron exchange
    o["ONEWAY"] = False
    if o.get("SOURCE") is None:
        print("No SOURCE file specified, creating a new one.")
        o["SOURCE"] = select_source_defaults(copy.deepcopy(options), option_descriptions)
    try:
        o["SOURCE_START"], note = resolve_source_start(
            o["SOURCE"], o.get("SOURCE_START"), get_mtime(options["inp"]["SOURCE"]))
    except ValueError as e:
        raise ValueError(f"\033[31m{e}{RESET}") from None
    if note:
        print(f"{YELLOW}{note}{RESET}")
    START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(o["start_time"], o["stop_time"])
    o["START_YEAR"] = START_YEAR
    o["START_DAY"] = START_DAY
    o = options.setdefault("job", {})
    hpc_platform = options["simulation"]["hpc_system"]
    for key in ("queue", "walltime", "job_priority", "group_list"):
        if key in option_descriptions["job"].get(hpc_platform, {"queue": 1, "walltime": 1}):
            if engage_options.get(key) is not None:
                o[key] = engage_options[key]
    on = options["job"].setdefault("resource", {})
    if hpc_platform == "aitken":
        on["model"] = engage_options["model"]
    select_default, ncpus_default, mpiprocs_default = select_resource_defaults(options, option_descriptions)
    for key, value in (("select", select_default), ("ncpus", ncpus_default), ("mpiprocs", mpiprocs_default)):
        if on.get(key) is None:
            on[key] = value
    if hpc_platform == "aitken":
        for key in ("local_modules", "other_job"):
            if options["job"].get(key) is None:
                options["job"][key] = option_descriptions["job"]["aitken"][key]["default"]
    o["nprocs"] = int(on["select"]) * int(on["mpiprocs"])
    if "project_code" in option_descriptions["job"].get(hpc_platform, {}):
        o["project_code"] = engage_options["project_code"]
        if o["project_code"] in (None, "", "None"):
            raise ValueError(f"\033[31mpbs.account_name (engage) is required on {hpc_platform}: it is the "
                             f"TIE-GCM jobs' project_code (#PBS -A).{RESET}")
    else:
        o.pop("project_code", None)
    note_engage_owned(saved, options)
    return options


# TIE-GCM options set by engage in a coupled run.
ENGAGE_OWNED_KEYS = (("simulation", "job_name"), ("simulation", "hpc_system"),
                     ("model", "specification", "horires"), ("inp", "ONEWAY"),
                     ("inp", "start_time"), ("inp", "stop_time"),
                     ("inp", "secondary_start_time"), ("inp", "secondary_stop_time"), ("inp", "segment"),
                     ("job", "queue"), ("job", "walltime"), ("job", "job_priority"), ("job", "group_list"),
                     ("job", "project_code"), ("job", "resource", "model"))


def _path_value(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def note_engage_owned(saved, options):
    """Print one NOTE for engage-owned keys whose saved value engage replaced; returns their names."""
    changed = []
    for path in ENGAGE_OWNED_KEYS:
        old, new = _path_value(saved, path), _path_value(options, path)
        if old in (None, "", "None") or str(old) == str(new):
            continue
        changed.append((".".join(path), old, new))
    if changed:
        print(f"{YELLOW}NOTE: engage owns these TIE-GCM options and replaced the -to JSON's values "
              f"(change them in the engage options): "
              + ", ".join(f"{name} {old!r} -> {new!r}" for name, old, new in changed) + f".{RESET}")
    return [name for name, _, _ in changed]


def get_engage_start_time(datetime_str, seconds):
    """`seconds` before `datetime_str`, rounded down to midnight."""
    dt = datetime.fromisoformat(datetime_str)

    new_dt = dt - timedelta(seconds=seconds)

    if new_dt.time() != datetime.min.time():
        new_dt = datetime.combine(new_dt.date(), datetime.min.time())

    return new_dt.isoformat()

def coupled_inp_defaults(inp, horires):
    """COUPLED_CAPS plus the `horires` GSWM files; GSWM is None where another lower-boundary
    source is set or the model refuses it."""
    out = dict(COUPLED_CAPS)
    other = lbc_other_set(inp)
    for key, pattern in GSWM_PATTERNS.items():
        allowed = gswm_allowed(key, horires) and not other
        out[key] = (find_file(pattern.format(horires=float(horires)), tiegcm_env("TIEGCMDATA"))
                    if allowed else None)
    return out


def apply_coupled_defaults(inp, horires, keep_unset=()):
    """Fill unset inp keys (except `keep_unset`) with the coupled defaults; returns the keys filled."""
    filled = []
    for key, value in coupled_inp_defaults(inp, horires).items():
        if key in keep_unset:
            continue
        if is_unset(inp.get(key)):
            inp[key] = value
            if value is not None:
                filled.append(key)
    return filled


def coupled_gswm(inp, horires_standalone, horires_coupled):
    """Switch set GSWM files to the coupled grid's; returns a NOTE for replaced non-defaults, or None."""
    if float(horires_coupled) == float(horires_standalone):
        return None
    spinup = coupled_inp_defaults(inp, horires_standalone)
    coupled = coupled_inp_defaults(inp, horires_coupled)
    replaced = []
    for key in GSWM_PATTERNS:
        if is_unset(inp.get(key)):
            continue
        if inp[key] != spinup[key]:
            replaced.append(f"{key} {inp[key]} -> {coupled[key]}")
        inp[key] = coupled[key]
    if not replaced:
        return None
    note = (f"NOTE: the coupled run is on the {float(horires_coupled)}-degree grid: "
            + "; ".join(replaced))
    print(f"{YELLOW}{note}{RESET}")
    return note


def _spinup_sechist(start_time, stop_time, preferred_sec_s, step_s):
    """Spin-up SECHIST in seconds that divides every segment, including a partial last one
    (TIE-GCM aborts otherwise)."""
    span = int((datetime.fromisoformat(stop_time) - datetime.fromisoformat(start_time)).total_seconds())
    sec = math.gcd(preferred_sec_s, span)
    if sec <= 0 or sec % step_s != 0:
        raise ValueError(
            f"coupled_start={stop_time!r} is too finely off-midnight: the spin-up SECHIST "
            f"cadence derived to divide the {span}s span is {sec}s, which is not a multiple "
            f"of STEP {step_s}s. Use a coupled start aligned to a coarser (whole-minute / "
            f"whole-hour) boundary.")
    return sec


VERIFIED_COUPLED_MACHINES = ("derecho", "aitken")


def experimental_machine_warning(hpc_system):
    if hpc_system in VERIFIED_COUPLED_MACHINES:
        return None
    return (f"WARNING: a coupled (engage GRT) TIE-GCM run on {hpc_system} is EXPERIMENTAL: it has never "
            f"run there, and the job scripts (from config/machines.yaml) are unchecked on that machine.")


def ensure_job_block(options, option_descriptions, engage):
    """Add a default job block where the machine has no job prompts (e.g. linux)."""
    job = options.get("job")
    if isinstance(job, dict) and job.get("resource"):
        return options
    job = options["job"] = dict(job or {})
    select, ncpus, mpiprocs = select_resource_defaults(options, option_descriptions)
    job["resource"] = {"select": select, "ncpus": ncpus, "mpiprocs": mpiprocs}
    job["nprocs"] = int(select) * int(mpiprocs)
    for key in ("queue", "walltime"):
        job.setdefault(key, (engage or {}).get(key))
    return options


class CoupledOptions(dict):
    """The coupled TIE-GCM options; .replay_options holds the input options to save for a -to
    replay (saving the coupled ones would make the replay regrid its own product)."""
    replay_options = None


def engage_run(options, debug, coupling, engage):
    """Write the TIE-GCM spin-up jobs and coupled namelists. Returns (CoupledOptions, spin-up
    .pbs files, coupled .inp files), paths relative to the run directory."""
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = apply_system_config(json.load(f), socket.gethostname())
    options = copy.deepcopy(options)
    horires_standalone = float(engage["horires"])
    unset = [k for k in engage.get("coupled_unset", ()) if is_unset(options["inp"].get(k))]
    apply_coupled_defaults(options["inp"], horires_standalone, unset)
    options_replay = copy.deepcopy(options)
    if unset:
        # keys answered 'none' stay unset on replay
        options_replay["coupled_unset"] = unset
    warning = experimental_machine_warning(options["simulation"]["hpc_system"])
    if warning:
        print(f"{YELLOW}{warning}{RESET}")
    ensure_job_block(options, option_descriptions, engage)
    options_standalone = copy.deepcopy(options)
    options_coupling = copy.deepcopy(options)
    options_standalone["simulation"]["job_name"] = f'{engage["job_name"]}-tiegcm-standalone'
    options_standalone["inp"]["stop_time"] = engage["coupled_start_time"]
    if float(engage['horires']) <= 2.5:
        standalone_segment_length_days = 7
    else:
        standalone_segment_length_days = 3
    options_standalone["inp"]["segment"] = f'{standalone_segment_length_days} 0 0 0'
    options_standalone["model"]["data"]["workdir"] = os.path.join(engage["parentdir"],"tiegcm_standalone")
    options_standalone["model"]["data"]["histdir"] = os.path.join(engage["parentdir"],"tiegcm_standalone")

    in_prim = options_standalone["inp"]["SOURCE"]
    out_prim = f'{options_standalone["model"]["data"]["workdir"]}/{options_standalone["simulation"]["job_name"]}_prim.nc'
    options_standalone["inp"]["SOURCE"] = out_prim
    # the spin-up uses its grid's default STEP, not dtCouple
    vertres_standalone, mres_standalone, nres_grid_standalone, STEP_standalone = resolution_solver(horires_standalone)
    options_standalone["inp"]["STEP"] = STEP_standalone
    # one primary history per segment, at its end, for the next segment to start from
    segment_s = standalone_segment_length_days * 86400
    sec_s = _spinup_sechist(engage["start_time"], engage["coupled_start_time"], 3600, int(STEP_standalone))
    options_standalone["inp"]["PRIHIST"] = options_standalone["inp"]["segment"]
    options_standalone["inp"]["MXHIST_PRIM"] = 1
    options_standalone["inp"]["SECHIST"] = " ".join(map(str, seconds_to_dhms(sec_s)))
    options_standalone["inp"]["MXHIST_SECH"] = segment_s // sec_s
    # plan both runs before writing anything, so a refused value leaves no partial output
    _, _, pristop_times, standalone_last_prim = segment_plan(options_standalone, options_standalone["simulation"]["job_name"])
    if not pristop_times:
        raise ValueError(
            f"The TIE-GCM spin-up window is empty: it starts at the midnight at or before the coupled "
            f"start minus gr_warm_up_time + gcm_spin_up_time ({engage['start_time']}) and ends at the "
            f"coupled start ({engage['coupled_start_time']}). Set coupling.gcm_spin_up_time > 0.")
    options_coupling["model"]["data"]["modelexe"] = options_coupling["model"]["data"]["coupled_modelexe"]
    horires_coupling = float(engage["horires_coupled"])
    options_coupling["model"]["specification"]["horires"] = horires_coupling
    vertres_coupling, mres_coupling, nres_grid_coupling, STEP_coupling = resolution_solver(horires_coupling,engage)
    options_coupling["model"]["specification"]["vertres"] = vertres_coupling
    options_coupling["model"]["specification"]["mres"] = mres_coupling
    options_coupling["model"]["specification"]["nres_grid"] = nres_grid_coupling
    options_coupling["inp"]["STEP"] = STEP_coupling
    if coupled_handoff_regrids(engage):
        # the last spin-up job regrids its final history onto the coupled grid into this file
        SOURCE_coupling = os.path.join(engage["parentdir"], f'{engage["job_name"]}_prim.nc')
    else:
        SOURCE_coupling = standalone_last_prim
    options_coupling["inp"]["SOURCE"] = SOURCE_coupling
    # SOURCE_START is the mtime nc_wrhist stores, which wraps past year end (Jan 1 00:00 -> 1 0 0 0)
    last_segment_start = segment_time(options_standalone["inp"]["start_time"], options_standalone["inp"]["stop_time"],
                                      parse_segment(options_standalone["inp"]["segment"]))[-1][0]
    last_segment_year = parse_run_datetime(last_segment_start, "spin-up segment start").year
    last_pristop_arr = list(map(int, pristop_times[-1].split()))
    options_coupling["inp"]["SOURCE_START"] = " ".join(map(str, nc_wrhist_mtime(last_pristop_arr, last_segment_year)))
    options_coupling["simulation"]["job_name"] = f'{engage["job_name"]}'
    options_coupling["inp"]["start_time"] = engage["coupled_start_time"]
    options_coupling["inp"]["secondary_start_time"] = engage["coupled_start_time"]
    options_coupling["inp"]["PRIHIST"] = " ".join(str(i) for i in engage["segment"])
    options_coupling["inp"]["MXHIST_PRIM"] = 1
    options_coupling["inp"]["SECHIST"] = " ".join(str(i) for i in seconds_to_dhms(engage["voltron_dtOut"]))
    options_coupling["inp"]["MXHIST_SECH"] = int(engage["segment_seconds"]/engage["voltron_dtOut"])
    options_coupling["inp"]["LABEL"] = f'{engage["job_name"]}_{horires_coupling}x{vertres_coupling}'
    options_coupling["inp"]["ONEWAY"] = False
    coupled_gswm(options_coupling["inp"], horires_standalone, horires_coupling)
    segment_plan(options_coupling, options_coupling["simulation"]["job_name"])
    # sized for the coupled grid, not the spin-up's job.resource
    coupled_resource = dict(zip(RESOURCE_KEYS, select_resource_defaults(options_coupling, option_descriptions,
                                                                        queue_sized=False)))
    for key in RESOURCE_KEYS:
        options_coupling["job"]["resource"][key] = coupled_resource[key]
    nprocs_coupling = int(coupled_resource["mpiprocs"])*int(coupled_resource["select"])
    options_coupling["job"]["nprocs"] = nprocs_coupling
    problems = job_resource_problems(options_standalone, where="job")
    error = ntask_error(nprocs_coupling, horires_coupling, nres_grid_coupling,
                        where=f"the coupled TIE-GCM ranks (select {coupled_resource['select']} x mpiprocs "
                              f"{coupled_resource['mpiprocs']})")
    if error:
        problems.append(error)
    problems += [p.replace("job.resource", "the coupled TIE-GCM ranks (job.resource)", 1)
                 if p.startswith("job.resource") else p
                 for p in validate_options(options_coupling, option_descriptions, generated=("inp.SOURCE",))
                 if p not in problems and not (error and p.startswith("job.nprocs"))]
    if problems:
        raise ValueError("\033[31mThe TIE-GCM jobs of this coupled run would be rejected:\n  - "
                         + "\n  - ".join(problems) + RESET)
    standalone_planned = segment_jobs(options_standalone, options_standalone["simulation"]["job_name"], True, engage)
    coupling_planned = segment_jobs(options_coupling, options_coupling["simulation"]["job_name"], False, engage)
    os.makedirs(options_standalone["model"]["data"]["workdir"], exist_ok=True)
    interpic (in_prim,float(horires_standalone),float(vertres_standalone),float(options_standalone['model']['specification']['zitop']),out_prim)
    standalone_inp_files,standalone_pbs_files, standalone_log_files,pristart_times, pristop_times, standalone_last_prim=write_segment_jobs(standalone_planned, options_standalone["simulation"]["job_name"])
    coupling_inp_files,coupling_pbs_files, coupling_log_files, pristart_times, pristop_times, _ = write_segment_jobs(coupling_planned, options_coupling["simulation"]["job_name"])

    options_coupling = CoupledOptions(options_coupling)
    options_coupling.replay_options = json_paths(options_replay, engage["parentdir"])
    run_relative = lambda p: os.path.relpath(p, engage["parentdir"]) if p else p      # noqa: E731
    return (options_coupling, [run_relative(p) for p in standalone_pbs_files],
            [run_relative(p) for p in coupling_inp_files])
