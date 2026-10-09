"""Replay of a saved options JSON (tiegcmrun -o) and the option checks every entry path runs.

Omitted keys are filled as the prompts would fill them; --rederive drops groups of derived keys
so they are filled again."""

import copy
import os
import re

import math

from misc import (resolution_solver, mres_to_nres_grid, find_file, GSWM_PATTERNS, gswm_allowed,
                  he_coefs_pattern, unknown_other_input_keys, other_input_error, secflds_list,
                  get_mtime, select_resource_defaults, select_source_defaults, history_dir,
                  default_make_file, MAKE_FRAGMENTS, tiegcm_env, coerce, is_unset)
from namelist_solver import (inp_pri_date, inp_prihist, inp_sechist, inp_mxhist, inp_pri_out,
                             inp_sec_out, inp_sec_date, parse_segment, forcing_missing,
                             forcing_conflicts, forcing_refusal, source_time_of_day_error,
                             resolve_source_start, source_start_hms_error, window_problems,
                             lbc_conflicts, lbc_other_set, cadence_problems)

RED = '\033[31m'
YELLOW = '\033[33m'
RESET = '\033[0m'

REDERIVE_GROUPS = {
    "job": "the job block, model.data.make and model.data.tgcmdata, from config/machines.yaml for "
           "simulation.hpc_system (a saved project_code / group_list is kept)",
    "dates": "START_YEAR, START_DAY, PRISTART, PRISTOP, the secondary-history window, the history "
             "cadences and file counts and OUTPUT / SECOUT, from start_time, stop_time and segment",
    "names": "LABEL, OUTPUT, SECOUT and model.data.log_file, from job_name and the grid",
    "resolution": "vertres, mres, nres_grid, STEP and the GSWM / HE_COEFS files, from horires",
    "source": "SOURCE and SOURCE_START: the nearest seasonal file at solar_flux_level and its "
              "history at the run start",
}

# Generated OUTPUT/SECOUT forms: group 1 is the history directory, group 2 the run name.
_GENERATED_OUTPUT = [
    re.compile(r"^'(.+)/([^/]+)_prim_\d+\.nc','to','\1/\2_prim_\d+\.nc','by','1'$"),
    re.compile(r"^'(.+)/([^/]+)_temp_\d+\.nc' , '\1/\2_prim_\d+\.nc'$"),
]
_GENERATED_SECOUT = [
    re.compile(r"^'(.+)/([^/]+)_sech_\d+\.nc','to','\1/\2_sech_\d+\.nc','by','1'$"),
    re.compile(r"^'(.+)/([^/]+)_sech_\d+\.nc'$"),
]
# <job_name>_<horires>x<vertres>
_RUN_NAME = re.compile(r"^(.+)_(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)$")
_GSWM_RES = re.compile(r"_(\d+(?:\.\d+)?)d_99km")

# Derived inp keys, in dependency order.
_TIME_KEYS = ("START_YEAR", "START_DAY", "PRISTART", "PRISTOP")
_DERIVED_KEYS = (("STEP",) + _TIME_KEYS
                 + ("secondary_start_time", "secondary_stop_time", "PRIHIST", "MXHIST_PRIM", "OUTPUT",
                    "SECHIST", "MXHIST_SECH", "SECSTART", "SECSTOP", "SECOUT", "LABEL",
                    "HE_COEFS_NCFILE") + tuple(GSWM_PATTERNS))
_REQUIRED_KEYS = ("SOURCE", "SOURCE_START")

_GROUP_INP_KEYS = {
    "dates": _TIME_KEYS + ("secondary_start_time", "secondary_stop_time", "PRIHIST", "MXHIST_PRIM",
                           "SECHIST", "MXHIST_SECH", "SECSTART", "SECSTOP", "OUTPUT", "SECOUT"),
    "names": ("LABEL", "OUTPUT", "SECOUT"),
    "resolution": ("STEP", "HE_COEFS_NCFILE") + tuple(GSWM_PATTERNS),
    "source": ("SOURCE", "SOURCE_START"),
}

ACCOUNT_KEYS = ("project_code", "group_list")


# Job keys an older tiegcmrun wrote that config/machines.yaml now owns.
OLD_JOB_KEYS = ("modules", "moduledir", "mpi_command", "account_name")


def old_format_problem(options):
    """An error message for keys of an older tiegcmrun that would be ignored, else None."""
    old = [f"{k} is no longer read" for k in ("provenance",) if k in options]
    job = options.get("job")
    if isinstance(job, dict):
        old += [f"job.{k} is no longer read (the coupled TIE-GCM resource is always derived from the "
                f"coupled grid)" for k in ("coupled_resource", "coupled_resource_for") if k in job]
        old += [f"job.{k} is no longer read ("
                f"{'the account is job.project_code' if k == 'account_name' else 'from config/machines.yaml'})"
                for k in OLD_JOB_KEYS if k in job]
    if not old:
        return None
    return (f"the options JSON was written by an older tiegcmrun: {'; '.join(old)}. Remove "
            f"{'it' if len(old) == 1 else 'them'} from the JSON, or generate a new one.")


def parse_rederive(value):
    """The --rederive groups of value ('dates,names', 'all', 'list') as a tuple."""
    groups = [g.strip().lower() for g in str(value).split(",") if g.strip()]
    if groups == ["list"]:
        return ("list",)
    if "all" in groups:
        return tuple(REDERIVE_GROUPS)
    unknown = [g for g in groups if g not in REDERIVE_GROUPS]
    if unknown or not groups:
        raise ValueError(f"unknown --rederive group {', '.join(unknown) or repr(value)} "
                         f"(groups: {', '.join(REDERIVE_GROUPS)}, all, list)")
    return tuple(dict.fromkeys(groups))


def rederive_groups_text():
    return "\n".join(f"{name:<11} {what}" for name, what in REDERIVE_GROUPS.items())


def _ints(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [int(v) for v in value]
    return [int(v) for v in str(value).replace(",", " ").split()]


def _dhms(values):
    return " ".join(str(v) for v in values)


def _set(value):
    return value is not None and str(value).strip().lower() not in ("", "none", "null")


def _run_name(options):
    spec = options["model"]["specification"]
    return f"{options['simulation']['job_name']}_{spec['horires']}x{spec['vertres']}"


def _generated_name(value, patterns):
    """(histdir, run_name) of a generated-form OUTPUT/SECOUT value, else None."""
    if not isinstance(value, str):
        return None
    for pattern in patterns:
        m = pattern.match(value.strip())
        if m:
            return m.group(1), m.group(2)
    return None


def rederive_pop(options, groups):
    """Remove the keys of each --rederive group in place; returns the removed job block's account."""
    groups = set(groups or ())
    inp = options.get("inp") if isinstance(options.get("inp"), dict) else {}
    model = options.get("model", {})
    data, spec = model.get("data", {}), model.get("specification", {})
    account = {}
    if "job" in groups:
        job = options.pop("job", None)
        if isinstance(job, dict):
            account = {k: job[k] for k in ACCOUNT_KEYS if _set(job.get(k))}
        data.pop("make", None)
        data.pop("tgcmdata", None)
    if "names" in groups:
        data.pop("log_file", None)
    if "resolution" in groups:
        for key in ("vertres", "mres", "nres_grid"):
            spec.pop(key, None)
    for group, keys in _GROUP_INP_KEYS.items():
        if group in groups:
            for key in keys:
                inp.pop(key, None)
    return account


def _derive(key, options, segment):
    """The value the prompts derive for inp[key] from the primary inputs."""
    inp = options["inp"]
    data = options["model"]["data"]
    horires = float(options["model"]["specification"]["horires"])
    run_name = _run_name(options)
    if key == "STEP":
        return resolution_solver(horires)[3]
    if key == "LABEL":
        return run_name
    if key in GSWM_PATTERNS:
        if not gswm_allowed(key, horires) or lbc_other_set(inp):
            return None
        return find_file(GSWM_PATTERNS[key].format(horires=horires), data["tgcmdata"])
    if key == "HE_COEFS_NCFILE":
        return find_file(he_coefs_pattern(horires), data["tgcmdata"])
    start, stop = inp.get("start_time"), inp.get("stop_time")
    if not (_set(start) and _set(stop)):
        raise ValueError(f"{RED}inp.{key} is missing from the options JSON and cannot be derived "
                         f"without inp.start_time and inp.stop_time: add inp.{key} (or the run "
                         f"window) to the JSON.{RESET}")
    START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(start, stop)
    if key in _TIME_KEYS:
        return {"START_YEAR": START_YEAR, "START_DAY": START_DAY,
                "PRISTART": _dhms(PRISTART), "PRISTOP": _dhms(PRISTOP)}[key]
    if key == "secondary_start_time":
        return start
    if key == "secondary_stop_time":
        return stop
    sec_start = inp["secondary_start_time"] if _set(inp.get("secondary_start_time")) else start
    sec_stop = inp["secondary_stop_time"] if _set(inp.get("secondary_stop_time")) else stop
    step = int(float(inp["STEP"])) if _set(inp.get("STEP")) else resolution_solver(horires)[3]
    if key == "PRIHIST":
        return _dhms(inp_prihist(_ints(inp["PRISTART"]), _ints(inp["PRISTOP"]), segment, step))
    if key == "MXHIST_PRIM":
        return inp_mxhist(start, stop, _ints(inp["PRIHIST"]), None, segment)[0]
    if key == "OUTPUT":
        return inp_pri_out(start, stop, _ints(inp["PRIHIST"]), inp["MXHIST_PRIM"], 0,
                           history_dir(data), run_name)[0]
    if key == "SECHIST":
        _, _, window_start, window_stop = inp_pri_date(sec_start, sec_stop)
        return _dhms(inp_sechist(window_start, window_stop, segment, step))
    if key == "MXHIST_SECH":
        return inp_mxhist(start, stop, _ints(inp["SECHIST"]), None, segment)[0]
    if key in ("SECSTART", "SECSTOP"):
        SECSTART, SECSTOP = inp_sec_date(sec_start, sec_stop, _ints(inp["SECHIST"]))
        return _dhms(SECSTART if key == "SECSTART" else SECSTOP)
    if key == "SECOUT":
        return inp_sec_out(sec_start, sec_stop, _ints(inp["SECHIST"]), inp["MXHIST_SECH"], 0,
                           history_dir(data), run_name)[0]
    raise KeyError(key)


def _pristart(inp):
    """PRISTART from start_time/stop_time, else the saved PRISTART (benchmark JSON)."""
    if _set(inp.get("start_time")) and _set(inp.get("stop_time")):
        return inp_pri_date(inp["start_time"], inp["stop_time"])[2]
    return _ints(inp["PRISTART"]) if _set(inp.get("PRISTART")) else None


def _derive_source(options, rederive=True):
    """The prompts' default SOURCE."""
    inp = options["inp"]
    what = "--rederive source:" if rederive else "The options JSON has no inp.SOURCE:"
    if not (_set(inp.get("start_time")) and _set(inp.get("solar_flux_level"))):
        need = "--rederive source needs" if rederive else "The options JSON has no inp.SOURCE: its default needs"
        raise ValueError(f"{RED}{need} inp.start_time and inp.solar_flux_level in the options "
                         f"JSON.{RESET}")
    data_dir = options["model"]["data"].get("tgcmdata") or tiegcm_env("TIEGCMDATA")
    source = select_source_defaults(options, None, data_dir)
    if source is None:
        raise ValueError(f"{RED}{what} no seasonal SOURCE file for this run under "
                         f"{data_dir}: set inp.SOURCE in the JSON.{RESET}")
    return source


def _derive_source_start(inp):
    """SOURCE's first history at the run start's time of day."""
    try:
        mtimes = get_mtime(inp["SOURCE"])
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as e:
        raise ValueError(f"{RED}inp.SOURCE_START is not in the options JSON and cannot be derived: "
                         f"SOURCE {inp.get('SOURCE')} cannot be read ({e}).{RESET}") from None
    return resolve_source_start(inp["SOURCE"], None, mtimes, pristart=_pristart(inp))[0]


def machine_data_dir(machine):
    """machines.yaml tgcmdata of machine, else $TIEGCMDATA when this host is it, else None."""
    from output_solver import CONFIG_DIR
    from jobgen import JobGen
    import socket
    m = _machines().get(machine, {})
    if m.get("tgcmdata"):
        return m["tgcmdata"]
    try:
        if JobGen(CONFIG_DIR, hostname=socket.gethostname()).name == machine:
            return tiegcm_env("TIEGCMDATA")
    except ValueError:
        pass
    return None


def _fill_omitted(options, option_descriptions, segment, groups=(), account=None):
    """Fill every key the JSON omits, so the template never renders an empty `KEY = `.
    Returns the names of the filled keys."""
    filled = []
    sim = options["simulation"]
    data = options["model"]["data"]
    spec = options["model"]["specification"]
    for section, name, keys in ((sim, "simulation", ("job_name", "hpc_system")),
                                (spec, "model.specification", ("horires", "zitop"))):
        descs = option_descriptions[name] if name == "simulation" else option_descriptions["model"]["specification"]
        for key in keys:
            if key not in section:
                section[key] = copy.deepcopy(descs[key].get("default"))
                filled.append(f"{name}.{key}")
    if "modeldir" not in data:
        data["modeldir"] = tiegcm_env("TIEGCMHOME")
        filled.append("model.data.modeldir")
    if not any(k in data for k in ("parentdir", "execdir", "workdir", "histdir")):
        data["parentdir"] = "."
        filled.append("model.data.parentdir")
    inp = options.get("inp")
    if isinstance(inp, dict) and inp:
        window = (_set(inp.get("start_time")) and _set(inp.get("stop_time"))) \
            or (_set(inp.get("PRISTART")) and _set(inp.get("PRISTOP")))
    else:
        window = _set(data.get("input_file"))
    if not window:
        raise ValueError(f"{RED}The options JSON has no value for inp.start_time, inp.stop_time: add "
                         f"them to the JSON.{RESET}")
    hpc = sim["hpc_system"]
    vertres, mres, _, _ = resolution_solver(float(spec["horires"]))
    for key, value in (("vertres", vertres), ("mres", mres)):
        if key not in spec:
            spec[key] = value
            filled.append(f"model.specification.{key}")
    if "nres_grid" not in spec:
        spec["nres_grid"] = mres_to_nres_grid(spec["mres"])
        filled.append("model.specification.nres_grid")
    if "tgcmdata" not in data:
        data["tgcmdata"] = machine_data_dir(hpc)
        if data["tgcmdata"] is None:
            data["tgcmdata"] = tiegcm_env("TIEGCMDATA")
            _data_dir_note(hpc, data["tgcmdata"])
        filled.append("model.data.tgcmdata")
    if "make" not in data and data.get("modeldir") and default_make_file(data["modeldir"], hpc):
        data["make"] = default_make_file(data["modeldir"], hpc)
        filled.append("model.data.make")
    parent = data.get("parentdir")
    if "execdir" not in data:
        data["execdir"] = os.path.join(parent, "exec") if _set(parent) else "."
        filled.append("model.data.execdir")
    for key, sub in (("workdir", "stdout"), ("histdir", "hist")):
        if key not in data:
            data[key] = os.path.join(parent, sub) if _set(parent) else data["execdir"]
            filled.append(f"model.data.{key}")
    if "modelexe" not in data and _set(data.get("execdir")):
        data["modelexe"] = os.path.join(data["execdir"], "tiegcm.exe")
        filled.append("model.data.modelexe")
    if "log_file" not in data and _set(data.get("workdir")):
        data["log_file"] = os.path.join(data["workdir"], f"{options['simulation']['job_name']}.out")
        filled.append("model.data.log_file")
    if options.get("job") is None and hpc != "linux":
        options["job"] = derive_job_block(options, option_descriptions, account)
        print(f"{YELLOW}NOTE: the options JSON has no job block: derived for {hpc} from "
              f"config/machines.yaml.{RESET}")
    elif isinstance(options.get("job"), dict) and hpc != "linux":
        filled += _fill_job_keys(options, option_descriptions)
    inp = options.get("inp")
    if not inp:
        return filled
    od = option_descriptions["inp"]
    # A benchmark JSON has no start_time/stop_time; their absence marks it as one.
    missing = [k for k in od if k not in inp and k not in ("start_time", "stop_time")]
    no_gpi = gpi_absent(inp) and _set(inp.get("start_time")) and _set(inp.get("stop_time"))
    if "solar_flux_level" in missing:
        inp["solar_flux_level"] = copy.deepcopy(od["solar_flux_level"].get("default"))
    if "SOURCE" in missing:
        inp["SOURCE"] = _derive_source(options, "source" in groups)
    for key in [k for k in _DERIVED_KEYS if k in missing]:
        inp[key] = _derive(key, options, segment)
    if "SOURCE_START" in missing:
        inp["SOURCE_START"] = _derive_source_start(inp)
    for key in [k for k in missing if k not in _DERIVED_KEYS + _REQUIRED_KEYS]:
        inp[key] = copy.deepcopy(od[key].get("default"))
    if no_gpi:
        default_gpi(inp, inp["start_time"], inp["stop_time"])
    return filled + [f"inp.{k}" for k in missing]


def _fill_job_keys(options, option_descriptions):
    """Fill the job keys a saved job block omits from the default block; returns their names."""
    job = options["job"]
    default = _default_job_block(options, option_descriptions, job=job) or {}
    filled = []
    res = job.setdefault("resource", {}) if "resource" in default else job.get("resource")
    for key, value in (default.get("resource") or {}).items():
        if key not in res:
            res[key] = value
            filled.append(f"job.resource.{key}")
    for key, value in default.items():
        if key in ("resource", "nprocs") or key in job:
            continue
        job[key] = value
        filled.append(f"job.{key}")
    if "nprocs" in default and "nprocs" not in job:
        try:
            job["nprocs"] = int(res["select"]) * int(res["mpiprocs"])
        except (KeyError, TypeError, ValueError):
            job["nprocs"] = default["nprocs"]
        filled.append("job.nprocs")
    return filled


def _time_hints(options, segment=None):
    """Saved time keys that differ from those derived from the run window."""
    inp = options["inp"]
    if not (_set(inp.get("start_time")) and _set(inp.get("stop_time"))):
        return []   # benchmark JSON
    START_YEAR, START_DAY, PRISTART, PRISTOP = inp_pri_date(inp["start_time"], inp["stop_time"])
    expected = {"START_YEAR": [START_YEAR], "START_DAY": [START_DAY], "PRISTART": PRISTART,
                "PRISTOP": PRISTOP}
    hints = [f"{k} {_dhms(_ints(inp[k]))} (from start_time and stop_time: {_dhms(v)})"
             for k, v in expected.items() if _set(inp.get(k)) and _ints(inp[k]) != v]
    if segment is None:
        keys = ("SECSTART", "SECSTOP") + (() if names_hint(options) else ("OUTPUT", "SECOUT"))
        for key in keys:
            if not _set(inp.get(key)):
                continue
            want = _derive(key, options, None)
            saved = _dhms(_ints(inp[key])) if key in ("SECSTART", "SECSTOP") else str(inp[key]).strip()
            if saved != str(want).strip():
                hints.append(f"{key} {saved} (derived: {want})")
    return hints


def _resolution_hints(options):
    """Saved grid values and GSWM files that are not horires' own."""
    spec = options["model"]["specification"]
    horires = float(spec["horires"])
    vertres, mres, nres_grid, STEP = resolution_solver(horires)
    hints = [f"{k} {spec[k]} (horires gives {v})"
             for k, v in (("vertres", vertres), ("mres", mres), ("nres_grid", nres_grid))
             if _set(spec.get(k)) and float(spec[k]) != float(v)]
    inp = options.get("inp") or {}
    if _set(inp.get("STEP")) and int(float(inp["STEP"])) > STEP:
        hints.append(f"STEP {inp['STEP']} (above horires' default {STEP})")
    for key in GSWM_PATTERNS:
        m = _GSWM_RES.search(os.path.basename(str(inp.get(key)))) if _set(inp.get(key)) else None
        if m and float(m.group(1)) != horires:
            hints.append(f"{key} is the {m.group(1)}-degree file")
    he = os.path.basename(str(inp.get("HE_COEFS_NCFILE"))) if _set(inp.get("HE_COEFS_NCFILE")) else ""
    want = he_coefs_pattern(horires).strip("*")
    if he.startswith("he_coefs_") and want not in he:
        hints.append(f"HE_COEFS_NCFILE is {he} (horires gives {want}.nc)")
    return hints


def names_hint(options, label_only=False):
    """A hint when LABEL (and OUTPUT / SECOUT unless label_only) name another run, else None."""
    inp = options.get("inp") or {}
    old = None if label_only else (_generated_name(inp.get("OUTPUT"), _GENERATED_OUTPUT)
                                   or _generated_name(inp.get("SECOUT"), _GENERATED_SECOUT))
    old_name = old[1] if old else None
    if old_name is None and isinstance(inp.get("LABEL"), str) and _RUN_NAME.match(inp["LABEL"]):
        old_name = inp["LABEL"]
    if old_name is None or old_name == _run_name(options):
        return None
    what = "inp.LABEL" if label_only else "inp.LABEL / OUTPUT / SECOUT"
    return (f"NOTE: {what} {'names' if label_only else 'name'} the run {old_name}, not {_run_name(options)} (from job_name and the "
            f"grid): pass --rederive names to follow it.")


def _machines():
    from output_solver import CONFIG_DIR
    from jobgen import JobGen
    return JobGen(CONFIG_DIR, machine="linux").machines


def _queue_names(machine, option_descriptions):
    """A machine's queue names from machines.yaml, else from options_description.json."""
    m = _machines().get(machine, {})
    names = list(m.get("queues") or [])
    if m.get("queue_default") and m["queue_default"] not in names:
        names.append(m["queue_default"])
    return names or list((option_descriptions.get("job", {}).get(machine, {}).get("queue") or {})
                         .get("valids") or [])


def _machine_evidence(options, option_descriptions):
    """Saved values that belong to a machine other than simulation.hpc_system."""
    hpc = options["simulation"]["hpc_system"]
    od_job = option_descriptions.get("job", {})
    machines = _machines()
    evidence = []
    data = options.get("model", {}).get("data", {})
    make = data.get("make")
    if (_set(make) and os.path.basename(make) in MAKE_FRAGMENTS.values()
            and os.path.basename(make) != MAKE_FRAGMENTS.get(hpc)):
        evidence.append(f"model.data.make is {os.path.basename(make)}")
    tgcmdata = data.get("tgcmdata")
    other_data = {m.get("tgcmdata") for name, m in machines.items() if name != hpc and m.get("tgcmdata")}
    if _set(tgcmdata) and tgcmdata in other_data and tgcmdata != machines.get(hpc, {}).get("tgcmdata"):
        evidence.append("model.data.tgcmdata is another machine's")
    job = options.get("job")
    if isinstance(job, dict) and hpc in od_job:
        mine = set(od_job[hpc]) | set(od_job.get("_common", {}))
        others = set().union(*[set(s) for m, s in od_job.items() if m not in ("_common", hpc)])
        foreign = sorted(k for k in job if k in others - mine)
        if foreign:
            evidence.append(f"job has {', '.join(foreign)}")
        queues = _queue_names(hpc, option_descriptions)
        if queues and _set(job.get("queue")) and job["queue"] not in queues:
            evidence.append(f"job.queue {job['queue']!r} is not one of {hpc}'s queues")
    return evidence


def machine_hint(options, option_descriptions):
    """A hint when the job block was built for another machine, else None."""
    evidence = _machine_evidence(options, option_descriptions)
    if not evidence:
        return None
    return (f"NOTE: simulation.hpc_system is {options['simulation']['hpc_system']} but "
            f"{'; '.join(evidence)}: pass --rederive job to re-derive the machine settings.")


def account_problem(options, option_descriptions, where="job.project_code"):
    """An error message when the machine requires project_code and the job block lacks it."""
    hpc = options["simulation"]["hpc_system"]
    if "project_code" not in option_descriptions.get("job", {}).get(hpc, {}):
        return None
    job = options.get("job")
    if not isinstance(job, dict) or _set(job.get("project_code")):
        return None
    return (f"{where} is required on {hpc} (the PBS account, #PBS -A; config/machines.yaml has no "
            f"shared default): add \"project_code\": \"<your project>\" to the job block.")


def _queue_problem(options, option_descriptions):
    """An error message for a job.queue no machine in machines.yaml has, else None."""
    job = options.get("job")
    hpc = options["simulation"]["hpc_system"]
    if not (isinstance(job, dict) and _set(job.get("queue")) and _queue_names(hpc, option_descriptions)):
        return None
    known = set()
    for machine in _machines():
        known.update(_queue_names(machine, option_descriptions))
    if job["queue"] in known:
        return None
    return (f"job.queue {job['queue']!r} is not a queue of any machine in config/machines.yaml "
            f"({hpc}: {' | '.join(_queue_names(hpc, option_descriptions))}).")


def _namelist_problems(options):
    """Unknown other_input keys and non-migrating GSWM files at 5 degrees."""
    inp = options["inp"]
    problems = []
    other = inp.get("other_input")
    lines = other if isinstance(other, (list, tuple)) else ([other] if _set(other) else [])
    unknown = unknown_other_input_keys(lines, options["model"]["data"].get("modeldir"))
    if unknown:
        problems.append("inp." + other_input_error(unknown))
    horires = float(options["model"]["specification"]["horires"])
    for key in GSWM_PATTERNS:
        if _set(inp.get(key)) and not gswm_allowed(key, horires):
            problems.append(f"inp.{key} is set: the model refuses non-migrating GSWM files at "
                            f"5 degrees (set it to null)")
    return problems


def secflds_replay(options, where="inp.SECFLDS"):
    """Check and normalise inp.SECFLDS in place (misc.secflds_list); ValueError for bad names."""
    inp = options.get("inp") or {}
    value = inp.get("SECFLDS")
    if value is None or value == [None]:
        return
    names, notes = secflds_list(value, options.get("model", {}).get("data", {}).get("modeldir"),
                                where=where, short=True)
    if names != [None] and names != list(value):
        inp["SECFLDS"] = names
        print(f"{YELLOW}NOTE: {where} normalised: {'; '.join(dict.fromkeys(notes))}.{RESET}")


def _forcing_problems(inp):
    """Missing and conflicting forcing values."""
    problems = []
    missing = forcing_missing(inp)
    if missing:
        problems.append(forcing_refusal(missing, "the missing values in the inp section of the options JSON"))
    problems += [f"inp.{why}" for _, why in forcing_conflicts(inp)]
    return problems


def _source_problems(inp):
    """SOURCE / SOURCE_START problems; skipped when SOURCE cannot be read."""
    source = inp.get("SOURCE")
    if not (_set(source) and os.path.isfile(str(source))):
        return []
    try:
        mtimes = get_mtime(source)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError):
        return []
    pristart = _pristart(inp)
    if pristart is not None:
        tod = source_time_of_day_error(source, mtimes, pristart, inp.get("start_time") if _set(inp.get("start_time")) else None)
        if tod:
            return [tod[1]]
    try:
        _, note = resolve_source_start(source, inp.get("SOURCE_START"), mtimes, pristart=pristart)
    except ValueError as e:
        return [str(e)]
    if note:
        return [f"inp.SOURCE_START {inp['SOURCE_START']} is not a history of SOURCE {source}: delete "
                f"inp.SOURCE_START (or pass --rederive source) to take SOURCE's history"]
    if pristart is not None:
        err = source_start_hms_error(inp["SOURCE_START"], pristart)
        if err:
            return [f"inp.{err}"]
    return []


def _default_job_block(options, option_descriptions, account=None, job=None):
    """The default job block of the machine sized for the queue, or None."""
    from output_solver import queue_walltime_default
    hpc = options["simulation"]["hpc_system"]
    od = option_descriptions.get("job", {}).get(hpc)
    if not od:
        return None
    given = job or {}
    out = {}
    for key, desc in od.items():
        if key == "resource":
            res = {k: copy.deepcopy(d.get("default")) for k, d in desc.items()}
            if _set((given.get("resource") or {}).get("model")):
                res["model"] = given["resource"]["model"]
            queue = given.get("queue") if _set(given.get("queue")) else (od.get("queue") or {}).get("default")
            probe = dict(options, job={"resource": res, "queue": queue})
            select, ncpus, mpiprocs = select_resource_defaults(probe, option_descriptions)
            res.update({k: v for k, v in (("select", select), ("ncpus", ncpus), ("mpiprocs", mpiprocs))
                        if k in res})
            out["resource"] = res
        elif key != "nprocs":
            out[key] = copy.deepcopy(desc.get("default"))
    queue = given.get("queue") if _set(given.get("queue")) else out.get("queue")
    if "walltime" in out and _set(queue):
        out["walltime"] = queue_walltime_default(hpc, queue) or out["walltime"]
    for key, value in (account or {}).items():
        if key in od:
            out[key] = value
    if "nprocs" in od:
        res = out.get("resource", {})
        out["nprocs"] = int(res["select"]) * int(res["mpiprocs"])
    return out


def derive_job_block(options, option_descriptions, account=None):
    """The prompts' default job block for the machine, keeping account; ValueError when a
    required setting has no default."""
    hpc = options["simulation"]["hpc_system"]
    job = _default_job_block(options, option_descriptions, account)
    if job is None:
        raise ValueError(f"The options JSON has no job block and simulation.hpc_system = {hpc} has "
                         f"no job settings in config/machines.yaml / options_description.json to "
                         f"derive one from: add a job block for {hpc} to the JSON.")
    od = option_descriptions["job"][hpc]
    unset = [f"job.{k}" for k in ("project_code", "queue", "walltime") if k in od and not _set(job.get(k))]
    if unset:
        raise ValueError(f"The options JSON has no job block, and {', '.join(unset)} "
                         f"{'has' if len(unset) == 1 else 'have'} no default for {hpc} in "
                         f"config/machines.yaml: add a job block with {', '.join(unset)} to the JSON.")
    return job


def ensure_sections(options):
    """Create missing top-level sections in place; ValueError for a section that is not an object."""
    for path in (("simulation",), ("model",), ("model", "data"), ("model", "specification")):
        node = options
        for key in path[:-1]:
            node = node[key]
        value = node.setdefault(path[-1], {})
        if not isinstance(value, dict):
            raise ValueError(f"{'.'.join(path)} must be a JSON object")
    return options


def read_options_json(path, shown=None):
    """The JSON object in path; ValueError naming the file otherwise."""
    import json
    shown = shown or path
    try:
        with open(path, "r", encoding="utf-8") as f:
            text = f.read()
    except FileNotFoundError:
        raise ValueError(f"{shown} does not exist") from None
    except UnicodeDecodeError:
        raise ValueError(f"{shown} is not valid JSON (not UTF-8 text)") from None
    except OSError as e:
        raise ValueError(f"{shown} cannot be read ({e.strerror or e})") from None
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"{shown} is not valid JSON (line {e.lineno} column {e.colno}: {e.msg})") from None
    if not isinstance(data, dict):
        raise ValueError(f"{shown} must be a JSON object")
    return data


# Sections that must be objects when present; inp and job may also be null.
_OBJECT_SECTIONS = (("simulation",), ("model",), ("model", "data"), ("model", "specification"))
_OBJECT_OR_NULL_SECTIONS = (("inp",), ("job",), ("job", "resource"))


def section_problem(options):
    """An error message for the first section that is present but not an object, else None."""
    for path, nullable in ([(p, False) for p in _OBJECT_SECTIONS]
                           + [(p, True) for p in _OBJECT_OR_NULL_SECTIONS]):
        node = options
        for key in path[:-1]:
            node = node.get(key) if isinstance(node, dict) else None
        if not isinstance(node, dict) or path[-1] not in node:
            continue
        value = node[path[-1]]
        if not isinstance(value, dict) and not (nullable and value is None):
            return f"{'.'.join(path)} must be a JSON object"
    return None


# Keys saved besides the options_description.json options.
_EXTRA_KEYS = ("simulation.node_type", "model.data.utildir")
_BOOKKEEPING_KEYS = ("coupled_unset",)
# Keys older tiegcmrun versions wrote, accepted so old files still replay.
_PREVIOUS_KEYS = ("inp.CADENCE", "inp.GSWM_data", "inp.start_date", "inp.stop_date",
                  "model.data.output_file", "job.module_list", "job.other", "job.coupled_resource",
                  "job.coupled_resource_for", "provenance", "pbs") + tuple(f"job.{k}" for k in OLD_JOB_KEYS)


def options_key_schema(option_descriptions):
    """The keys an options JSON may hold, as {key: sub-schema, or None for a leaf}."""
    od = option_descriptions
    job, resource = {}, {}
    for block in od.get("job", {}).values():
        for key, desc in (block or {}).items():
            job[key] = None
            if key == "resource" and isinstance(desc, dict):
                resource.update(dict.fromkeys(desc))
    job["resource"] = resource
    schema = {"simulation": dict.fromkeys(od["simulation"]),
              "model": {"data": dict.fromkeys(od["model"]["data"]),
                        "specification": dict.fromkeys(od["model"]["specification"])},
              "inp": dict.fromkeys(od["inp"]), "job": job}
    schema.update(dict.fromkeys(_BOOKKEEPING_KEYS))
    for path in _EXTRA_KEYS + _PREVIOUS_KEYS:
        *parents, leaf = path.split(".")
        node = schema
        for key in parents:
            node = node.setdefault(key, {})
        node.setdefault(leaf, None)
    return schema


def unknown_keys_problem(options, schema):
    """An error message naming keys not in schema, with a close match; keys starting with '_'
    are allowed."""
    import difflib

    def paths(sch, prefix=""):
        out = []
        for k, v in sch.items():
            out.append(prefix + k)
            if isinstance(v, dict):
                out += paths(v, f"{prefix}{k}.")
        return out

    everywhere = paths(schema)
    items = []

    def walk(node, sch, prefix):
        for k, v in node.items():
            if str(k).startswith("_"):
                continue
            if k not in sch:
                close = difflib.get_close_matches(str(k), list(sch), n=1)
                hint = (prefix + close[0]) if close else next(
                    (p for p in everywhere if p.rsplit(".", 1)[-1] == k), None)
                items.append(f"{prefix}{k}" + (f" (did you mean {hint}?)" if hint else ""))
            elif isinstance(sch[k], dict) and isinstance(v, dict):
                walk(v, sch[k], f"{prefix}{k}.")

    walk(options, schema, "")
    if not items:
        return None
    return f"unknown key{'' if len(items) == 1 else 's'} " + "; ".join(items)


_MANUAL_FORCING = ("KP", "POWER", "CTPOTEN", "F107", "F107A")


def gpi_absent(inp):
    """True when inp has no GPI_NCFILE key and no manual forcing value."""
    return "GPI_NCFILE" not in inp and not any(_set(inp.get(k)) for k in _MANUAL_FORCING)


def default_gpi(inp, start, stop):
    """Set inp.GPI_NCFILE to the newest bundled GPI file if it covers the window, else 'gen'."""
    from misc import select_latest_gpi, file_covers
    latest = select_latest_gpi(tiegcm_env("TIEGCMDATA"))
    inp["GPI_NCFILE"] = latest if latest and file_covers(latest, start, stop) else "gen"
    return inp["GPI_NCFILE"]


def fill_note(names):
    """Print one NOTE naming the filled keys."""
    if len(names) <= 8:
        text = f"{', '.join(names)}: using the default{'' if len(names) == 1 else 's'}."
    else:
        text = (f"{len(names)} keys ({', '.join(names[:3])}, ...): using the defaults; the full set is saved "
                f"in tiegcmrun_parameters.json.")
    print(f"{YELLOW}NOTE: the tiegcmrun options JSON omits {text}{RESET}")


def _data_dir_note(hpc, data_dir):
    print(f"{YELLOW}NOTE: {hpc} has no tgcmdata in config/machines.yaml: model.data.tgcmdata is "
          f"$TIEGCMDATA ({data_dir}); set it to {hpc}'s TIE-GCM data directory.{RESET}")


def fill_engage_omitted(options, option_descriptions, engage_owned=(), engage=None):
    """Fill omitted keys for engage's -to replay, leaving the keys engage owns or derives itself.
    Returns the filled keys."""
    inp = options.setdefault("inp", {})
    od = option_descriptions["inp"]
    derived = ("LABEL", "HE_COEFS_NCFILE")
    skip = (set(engage_owned) | set(_DERIVED_KEYS) | set(_REQUIRED_KEYS)
            | {"start_time", "stop_time", "segment"}) - set(derived)
    filled = []
    if engage is not None:
        data, spec = options["model"]["data"], options["model"]["specification"]
        parent = engage["parentdir"]
        for key, value in (("modeldir", lambda: tiegcm_env("TIEGCMHOME")),
                           ("modelexe", lambda: os.path.join(parent, "tiegcm.exe")),
                           ("coupled_modelexe", lambda: os.path.join(parent, "tiegcm.x"))):
            if key not in data:
                data[key] = value()
                filled.append(f"model.data.{key}")
        if "zitop" not in spec:
            spec["zitop"] = copy.deepcopy(option_descriptions["model"]["specification"]["zitop"].get("default"))
            filled.append("model.specification.zitop")
        if not data.get("tgcmdata"):
            data["tgcmdata"] = machine_data_dir(engage["hpc_system"])
            if data["tgcmdata"] is None:
                data["tgcmdata"] = tiegcm_env("TIEGCMDATA")
                _data_dir_note(engage["hpc_system"], data["tgcmdata"])
            filled.append("model.data.tgcmdata")
    no_gpi = engage is not None and gpi_absent(inp)
    for key in od:
        if key in inp or key in skip:
            continue
        inp[key] = _derive(key, options, None) if key in derived else copy.deepcopy(od[key].get("default"))
        filled.append(f"inp.{key}")
    if no_gpi:
        default_gpi(inp, engage["start_time"], engage["stop_time"])
    if filled:
        fill_note(filled)
    return filled


def _hints(options, option_descriptions):
    """One-line hints for derived values that no longer follow their primary inputs."""
    hints = []
    inp = options.get("inp")
    if inp:
        try:
            segment = parse_segment(inp.get("segment"))
        except ValueError:
            segment = None
        stale = _time_hints(options, segment)
        if stale:
            hints.append(f"NOTE: {', '.join(stale)}: pass --rederive dates to follow the run window.")
        hint = names_hint(options)
        if hint:
            hints.append(hint)
    stale = _resolution_hints(options)
    if stale:
        hints.append(f"NOTE: model.specification.horires is {options['model']['specification']['horires']} "
                     f"but {', '.join(stale)}: pass --rederive resolution to follow it.")
    if options.get("job") is not None and options["simulation"]["hpc_system"] != "linux":
        hint = machine_hint(options, option_descriptions)
        if hint:
            hints.append(hint)
    return hints


def prepare_replay(options, option_descriptions, rederive=()):
    """A copy of options ready to replay: --rederive groups removed, omitted keys filled, hints
    printed."""
    options = copy.deepcopy(options)
    ensure_sections(options)
    account = rederive_pop(options, rederive)
    inp = options.get("inp")
    try:
        segment = parse_segment(inp.get("segment")) if inp else None
    except ValueError:
        segment = None                        # refused by validate_options
    filled = _fill_omitted(options, option_descriptions, segment, rederive, account)
    if filled:
        fill_note(filled)
    try:
        hints = _hints(options, option_descriptions)
    except (ValueError, TypeError, KeyError):
        hints = []                            # a malformed value: validate_options names it
    for hint in hints:
        print(f"{YELLOW}{hint}{RESET}")
    return options


def settle_select(options):
    """Set job.resource.select to ceil(nprocs / mpiprocs) when they differ; returns the NOTE or None."""
    job = options.get("job")
    if not isinstance(job, dict):
        return None
    res = job.get("resource") or {}
    try:
        nprocs, mpiprocs = int(job["nprocs"]), int(res.get("mpiprocs", res.get("ncpus")))
        select = int(res["select"])
    except (KeyError, TypeError, ValueError):
        return None
    nodes = math.ceil(nprocs / mpiprocs)
    if nodes == select:
        return None
    res["select"] = nodes
    note = (f"NOTE: job.resource.select {select} -> {nodes} (nprocs {nprocs} / mpiprocs "
            f"{mpiprocs})")
    print(f"{YELLOW}{note}{RESET}")
    return note


def _active_job_descriptions(options, option_descriptions):
    od_job = option_descriptions.get("job", {})
    return {**(od_job.get("_common") or {}), **(od_job.get(options["simulation"].get("hpc_system")) or {})}


def _coerce_walk(options, option_descriptions, run_dir, generated, check_exists=True):
    """misc.coerce every set option value in place; returns the errors."""
    problems = []
    od = option_descriptions
    sections = [("simulation", options.get("simulation"), od.get("simulation", {})),
                ("model.data", (options.get("model") or {}).get("data"), od["model"]["data"]),
                ("model.specification", (options.get("model") or {}).get("specification"),
                 od["model"]["specification"]),
                ("inp", options.get("inp"), od.get("inp", {}))]
    job = options.get("job")
    if isinstance(job, dict):
        jd = _active_job_descriptions(options, od)
        sections.append(("job", job, jd))
        if isinstance(job.get("resource"), dict):
            sections.append(("job.resource", job["resource"], jd.get("resource") or {}))
    for prefix, values, descs in sections:
        if not isinstance(values, dict):
            continue
        for key, value in list(values.items()):
            desc = descs.get(key)
            name = f"{prefix}.{key}"
            if not isinstance(desc, dict) or "type" not in desc or is_unset(value):
                continue
            if name in generated or name == "inp.SECFLDS":
                continue                      # written by this run / checked by secflds_replay
            if name == "model.data.input_file" and options.get("inp"):
                continue                      # replaced by the generated deck
            if name == "job.queue":
                desc = {k: v for k, v in desc.items() if k != "valids"}   # valids are checked by _queue_problem
            try:
                values[key] = coerce(name, value, desc, run_dir, check_exists)
            except ValueError as e:
                problems.append(str(e))
    return problems


def other_cluster(machine):
    """True when machine is a cluster other than this host, so its files cannot be checked."""
    from output_solver import CONFIG_DIR
    from jobgen import JobGen
    import socket
    if machine in (None, "", "linux"):
        return False
    try:
        return JobGen(CONFIG_DIR, hostname=socket.gethostname()).name != machine
    except ValueError:
        return False


def validate_options(options, option_descriptions, run_dir=None, generated=()):
    """All problems with a run's options, checked before anything is written.

    Canonical values are written back into options; `generated` keys are skipped because this run
    writes or derives them itself."""
    from output_solver import job_resource_problems
    hpc = (options.get("simulation") or {}).get("hpc_system")
    problems = _coerce_walk(options, option_descriptions, run_dir, set(generated),
                            check_exists=not other_cluster(hpc))
    typed_ok = not problems

    def check(fn, *args):
        # A rule that cannot read a value the type check already refused stays quiet.
        try:
            return fn(*args) or []
        except (ValueError, TypeError, KeyError, IndexError, AttributeError) as e:
            return [str(e)] if typed_ok else []
    inp = options.get("inp")
    spec = (options.get("model") or {}).get("specification") or {}
    if isinstance(inp, dict) and inp:
        problems += check(_namelist_problems, options)
        try:
            secflds_replay(options)
        except ValueError as e:
            problems.append(str(e))
        problems += check(_forcing_problems, inp)
        if "inp.SOURCE" not in generated:
            problems += check(_source_problems, inp)
        problems += check(window_problems, inp)
        problems += check(lambda: [why for _, why in lbc_conflicts(inp, spec.get("horires"))])

        def cadence():
            segment = parse_segment(inp.get("segment"))
            step = inp.get("STEP") if _set(inp.get("STEP")) else resolution_solver(spec["horires"])[3]
            return cadence_problems({k: v for k, v in inp.items() if f"inp.{k}" not in generated}, segment, step)
        problems += check(cadence)
    if isinstance(options.get("job"), dict) and hpc != "linux":
        problems += [p for p in (account_problem(options, option_descriptions),
                                 _queue_problem(options, option_descriptions)) if p]
        if typed_ok:
            settle_select(options)
        problems += check(job_resource_problems, options)
    return list(dict.fromkeys(problems))
