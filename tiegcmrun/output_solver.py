"""Writes TIE-GCM namelists and job scripts, including segmented runs and the engage hand-off."""


import os
import re
import copy
import math
import shlex
from jinja2 import Template

from misc import (segment_time, resolution_solver, dhms_to_seconds, as_list, NAMELIST_ARRAYS, namelist_array,
                  job_path, history_dir, job_inp)
from namelist_solver import (inp_pri_date, inp_pri_out, inp_sec_date, inp_sec_out, parse_segment,
                             parse_run_datetime, fit_cadence, segment_mxhist, nc_wrhist_mtime,
                             RUN_DATETIME_FORMAT)
from jobgen import JobGen, account_value, user_projects


SUPPORT_FILES_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
OPTION_DESCRIPTIONS_FILE = os.path.join(SUPPORT_FILES_DIRECTORY, "options_description.json")
INP_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.inp")

CONFIG_DIR = os.path.join(SUPPORT_FILES_DIRECTORY, "config")


def _generic_job_descriptions(machine):
    """Job prompt descriptions for a machines.yaml machine options_description.json lacks."""
    def field(level, prompt, kind, default=None, **rule):
        return {"LEVEL": level, "type": kind, **rule, "prompt": prompt, "default": default,
                "description": None, "warning": None}
    return {
        "project_code": field("BENCH", "Project Code", "name"),
        "queue": field("INTERMEDIATE", f"{machine} queue name", "str"),
        "resource": {
            "select": field("EXPERT", "Number of nodes to request", "int", min=1),
            "ncpus": field("EXPERT", "Number of cores per node", "int", min=1),
            "mpiprocs": field("EXPERT", "Number of MPI ranks per node", "int", min=1),
        },
        "walltime": field("BASIC", "Requested wall time for each job segment (HH:MM:SS)", "hms"),
        "nprocs": field("EXPERT", "nprocs (nnodes * mpiprocs)", "int", min=1),
        "other_job": field("EXPERT", "Additional settings or env variables", "lines", [None]),
    }


def node_models(machine):
    """{PBS :model= tag: node type name} of a machines.yaml machine."""
    jg = JobGen(CONFIG_DIR, machine=machine)
    return {nt["model"]: name for name, nt in (jg.m.get("node_types") or {}).items() if nt.get("model")}


def apply_system_config(option_descriptions, hostname=None):
    """Seed the machine and job prompt defaults and valids from config/machines.yaml, in place."""
    jg = JobGen(CONFIG_DIR, hostname=hostname)
    hpc = option_descriptions["simulation"]["hpc_system"]
    hpc["valids"] = list(jg.machines)
    hpc["default"] = jg.name
    jobs = option_descriptions.setdefault("job", {})
    for name, m in jg.machines.items():
        if m.get("scheduler", "none") == "none":
            continue
        job = jobs.setdefault(name, _generic_job_descriptions(name))
        for key, yaml_key in (("project_code", "account_default"), ("queue", "queue_default"),
                              ("walltime", "walltime_default")):
            if key in job and m.get(yaml_key) is not None:
                job[key]["default"] = account_value(m[yaml_key]) if key == "project_code" else m[yaml_key]
                if job[key].get("valids") and m[yaml_key] not in job[key]["valids"]:
                    job[key]["valids"].append(m[yaml_key])
        if "project_code" in job and m.get("account_groups") and m.get("account_default") is None:
            choices = user_projects()
            if choices:
                job["project_code"]["choices"] = choices
        if "local_modules" in job and m.get("local_prefix"):
            job["local_modules"]["default"] = local_prefix_lines(m["local_prefix"], m.get("shell", "bash"))
        if "job_priority" in job and m.get("job_priorities"):
            job["job_priority"]["valids"] = list(m["job_priorities"])
        model = (job.get("resource") or {}).get("model")
        if model is not None and m.get("node_types"):
            models = node_models(name)
            model["valids"] = list(models)
            if model.get("default") not in models:
                default_nt = m["node_types"].get(m.get("node_type_default")) or next(iter(m["node_types"].values()))
                model["default"] = default_nt["model"]
    return option_descriptions


def local_prefix_lines(prefix, shell="tcsh"):
    """Job lines adding a machines.yaml local_prefix to the library, include and bin paths.

    tcsh aborts on an unset ${VAR}, even in a one-line if, hence the multi-line if/else."""
    if shell == "tcsh":
        L = [f"setenv PREFIX {prefix}"]
        for var, sub in (("LIBRARY_PATH", "lib"), ("LD_LIBRARY_PATH", "lib"), ("CPATH", "include")):
            L += [f"if ($?{var}) then", f"    setenv {var} $PREFIX/{sub}:${{{var}}}", "else",
                  f"    setenv {var} $PREFIX/{sub}", "endif"]
        return L + ["setenv PATH ${PATH}:$PREFIX/bin"]
    L = [f'export PREFIX="{prefix}"']
    for var, sub in (("LIBRARY_PATH", "lib"), ("LD_LIBRARY_PATH", "lib"), ("CPATH", "include")):
        L += [f'export {var}="$PREFIX/{sub}${{{var}:+:${var}}}"']
    return L + ['export PATH="${PATH}:$PREFIX/bin"']


def _hms_seconds(value):
    """'HH:MM:SS' -> seconds, or None."""
    m = re.fullmatch(r"\s*(\d+):(\d{1,2}):(\d{1,2})\s*", str(value))
    return int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3)) if m else None


def queue_limits(machine, queue):
    """{ncpus_max, walltime_max} of queue from machines.yaml, or {}."""
    queues = JobGen(CONFIG_DIR, machine=machine).m.get("queues") or {}
    if not isinstance(queues, dict):
        return {}
    return dict(queues.get(queue) or {})


def queue_walltime_default(machine, queue):
    """The machine's walltime_default, capped at the queue's walltime_max."""
    m = JobGen(CONFIG_DIR, machine=machine).m
    default, limit = m.get("walltime_default"), queue_limits(machine, queue).get("walltime_max")
    if default is None or (limit is not None and _hms_seconds(limit) < _hms_seconds(default)):
        return limit
    return default


def _node_cores(machine, res, node_type=None):
    """Cores per node of the job's node type."""
    if node_type is None and res.get("model") not in (None, "", "None"):
        node_type = node_models(machine).get(res["model"])
    return int(JobGen(CONFIG_DIR, machine=machine, node_type=node_type)._cores())


def job_resource_problems(options, where="job"):
    """Error messages for job resources the node, the queue limits or TIE-GCM's decomposition
    would reject."""
    from misc import ntask_error
    job = options.get("job") or {}
    machine = options["simulation"]["hpc_system"]
    problems = []
    res = job.get("resource") or {}
    limits = queue_limits(machine, job.get("queue")) if job.get("queue") else {}
    try:
        cores = _node_cores(machine, res, options["simulation"].get("node_type"))
    except (TypeError, ValueError):
        cores = None
    try:
        ncpus = res.get("ncpus")
        ncpus = int(ncpus) if ncpus not in (None, "", "None") else cores
        mpiprocs = res.get("mpiprocs", ncpus)
        mpiprocs = int(mpiprocs) if mpiprocs not in (None, "", "None") else ncpus
        nprocs = job.get("nprocs")
        if nprocs not in (None, "", "None"):
            select = math.ceil(int(nprocs) / mpiprocs)
        else:
            select = int(res.get("select")) if res.get("select") not in (None, "", "None") else None
    except (TypeError, ValueError, ZeroDivisionError):
        select = ncpus = mpiprocs = None
    if cores is not None and ncpus is not None and ncpus > cores:
        problems.append(f"{where}.resource.ncpus {ncpus} is more than the {cores} cores of a "
                        f"{machine} node.")
    if ncpus is not None and mpiprocs is not None and mpiprocs > ncpus:
        problems.append(f"{where}.resource.mpiprocs {mpiprocs} is more than ncpus {ncpus}.")
    if limits.get("ncpus_max") is not None and select is not None and select * ncpus > int(limits["ncpus_max"]):
        problems.append(f"{where}.resource: {select} nodes x {ncpus} cores = {select * ncpus} cores is more than "
                        f"queue {job['queue']!r} on {machine} allows ({limits['ncpus_max']}, "
                        f"config/machines.yaml queues): use fewer nodes or another queue.")
    wall, wall_max = _hms_seconds(job.get("walltime")), _hms_seconds(limits.get("walltime_max"))
    if wall is not None and wall_max is not None and wall > wall_max:
        problems.append(f"{where}.walltime {job['walltime']} is more than queue {job['queue']!r} on {machine} "
                        f"allows ({limits['walltime_max']}, config/machines.yaml queues): shorten it (and "
                        f"inp.segment) or use another queue.")
    nprocs = job.get("nprocs")
    if nprocs not in (None, "", "None"):
        spec = options["model"]["specification"]
        error = ntask_error(nprocs, spec["horires"], spec.get("nres_grid"), where=f"{where}.nprocs")
        if error:
            problems.append(error)
    return problems


def linux_run_commands(options, inp_files, log_files, relative=False):
    """mpirun commands for the namelists on a machine without a scheduler; they must run in the
    workdir. relative=True spells paths relative to it, else absolute."""
    jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"])
    nprocs = (options.get("job") or {}).get("nprocs") or jg._cores()
    groups = [{"mpi_tasks": int(nprocs)}]
    data = options["model"]["data"]
    name = (lambda p: job_path(p, data)) if relative else os.path.abspath   # noqa: E731
    return [jg.run_command(groups, name(data["modelexe"]), name(inp), name(log))
            for inp, log in zip(inp_files, log_files)]


# The OUTPUT forms inp_pri_out writes; the last file named is the last one the run fills.
_GENERATED_OUTPUT_FORMS = (re.compile(r"^\s*'([^']+)'\s*$"),
                           re.compile(r"^\s*'[^']+'\s*,\s*'([^']+)'\s*$"),
                           re.compile(r"^\s*'[^']+'\s*,\s*'to'\s*,\s*'([^']+)'\s*,\s*'by'\s*,\s*'1'\s*$"))


def last_primary_history(inp):
    """The last primary history file of a generated OUTPUT, or None for a hand-written one."""
    if not isinstance(inp, dict) or not inp.get("OUTPUT"):
        return None
    for form in _GENERATED_OUTPUT_FORMS:
        m = form.match(str(inp["OUTPUT"]))
        if m:
            return m.group(1).strip()
    return None


def completion_guard_lines(shell, log_file, stamp, last_prim=None):
    """Job lines after the model launch that exit non-zero unless TIE-GCM completed.

    The exit status alone is unreliable (shutdown() passes an uninitialised errorcode to
    mpi_abort), so the log and the age of the last primary history are checked too."""
    if shell == "tcsh":
        L = ["set rc = $status",
             "# Stop here (no hand-off, no next segment) unless tiegcm completed",
             "if ($rc != 0) then",
             '    echo "tiegcm failed (exit status $rc): stopping the job"',
             "    exit $rc",
             "endif",
             f'if ( {{ grep -q "MODEL SHUTDOWN" {log_file} }} ) then',
             f'    echo "MODEL SHUTDOWN in {log_file}: stopping the job"',
             "    exit 1",
             "endif"]
        if last_prim:
            L += [f"if ( ! -e {last_prim} ) then",
                  f'    echo "tiegcm did not write {last_prim}: stopping the job"',
                  "    exit 1",
                  "endif",
                  f'if ( "`find {last_prim} -newer {stamp}`" == "" ) then',
                  f'    echo "{last_prim} is older than this job: stopping the job"',
                  "    exit 1",
                  "endif"]
        return L
    L = ["rc=$?",
         "# Stop here (no hand-off, no next segment) unless tiegcm completed",
         'if [ $rc -ne 0 ]; then echo "tiegcm failed (exit status $rc): stopping the job"; exit $rc; fi',
         f'if grep -q "MODEL SHUTDOWN" {log_file}; then echo "MODEL SHUTDOWN in {log_file}: stopping the job"; exit 1; fi']
    if last_prim:
        L += [f'if [ ! -e {last_prim} ]; then echo "tiegcm did not write {last_prim}: stopping the job"; exit 1; fi',
              f'if [ -z "$(find {last_prim} -newer {stamp})" ]; then echo "{last_prim} is older than this job: stopping the job"; exit 1; fi']
    return L


def render_pbs(options):
    """Render a tiegcm job script from config/{machines,batch}.yaml and options.

    The job cds to its workdir by absolute path and names other files relative to it."""
    job = options["job"]
    data = options["model"]["data"]
    res = job.get("resource", {}) or {}
    node_type = options["simulation"].get("node_type")
    if node_type is None and res.get("model"):
        models = node_models(options["simulation"]["hpc_system"])
        if res["model"] not in models:
            raise ValueError(f"job.resource.model {res['model']!r} is not a node type of "
                             f"{options['simulation']['hpc_system']} in config/machines.yaml "
                             f"(known: {', '.join(models) or 'none'}).")
        node_type = models[res["model"]]
    jg = JobGen(CONFIG_DIR, machine=options["simulation"]["hpc_system"],
                node_type=node_type)
    if res.get("ncpus") not in (None, "", "None") and int(res["ncpus"]) != jg._cores():
        ncpus = int(res["ncpus"])
        jg._cores = lambda: ncpus
    task_groups = [{"mpi_tasks": int(job["nprocs"]),
                    "tasks_per_node": int(res.get("mpiprocs", res.get("ncpus", jg._cores()))),
                    "threads": int(res.get("ompthreads", 1))}]
    L = [jg.shebang(), "", "# Batch directives"]
    L += jg.directive_lines({
        "job_name": options["simulation"]["job_name"],
        "account": job.get("project_code"),
        "group_list": job.get("group_list"),
        "queue": job.get("queue"),
        "join_out": True,
        "priority": job.get("job_priority"),
    })
    L += jg.resource_lines(task_groups, job["walltime"], queue=job.get("queue"))
    L += ["", "# Set environment variables"]
    env = {"TGCMDATA": data["tgcmdata"], "TIEGCMDATA": data["tgcmdata"],
           "TIEGCMHOME": data["modeldir"]}
    env.update(jg.m.get("env", {}) or {})
    L += jg.env_lines(env)
    L += [ln for ln in as_list(job.get("other_job")) if ln is not None]
    L += ["", "# Load modules"]
    L += jg.module_lines()
    L += [ln for ln in as_list(job.get("local_modules")) if ln is not None]
    L += ["", "# Execute tiegcm in its work directory (the paths below are relative to it)",
          f"cd {shlex.quote(os.path.abspath(data['workdir']))}"]
    log = job_path(data["log_file"], data)
    stamp = f"{log}.start"
    L += [f"touch {stamp}"]
    L += [jg.run_command(task_groups, job_path(data["modelexe"], data),
                         job_path(data["input_file"], data), log)]
    last_prim = last_primary_history(options.get("inp"))
    L += completion_guard_lines(jg.m["shell"], log, stamp, last_prim)
    L += ["", "# Job chaining"]
    L += [ln for ln in as_list(job.get("job_chain")) if ln is not None]
    return "\n".join(L) + "\n"

def run_file(workdir, run_name, segment_number, ext):
    """<workdir>/<run_name>.<ext>, or <run_name>-NN.<ext> for segment NN (counted from 1)."""
    if segment_number is None:
        return os.path.join(workdir, f"{run_name}.{ext}")
    return os.path.join(workdir, f"{run_name}-{'{:02d}'.format(segment_number+1)}.{ext}")


def create_pbs_scripts(options, run_name, segment_number, content=None):
    """Write the job script (rendered here unless content is given) and return its path."""
    pbs_content = render_pbs(options) if content is None else content
    pbs_script = run_file(options["model"]["data"]["workdir"], run_name, segment_number, "pbs")
    with open(pbs_script, "w", encoding="utf-8") as f:
        f.write(pbs_content)
    return pbs_script

def create_inp_scripts(options, run_name, segment_number):
    """Write the namelist from template.inp and return its path."""
    with open(INP_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    # gfortran hits end-of-file on a deck whose last line is a bare '/' without a newline.
    template = Template(template_content, keep_trailing_newline=True)
    opt = copy.deepcopy(options)
    inp = opt.get("inp")
    if isinstance(inp, dict):
        for key in NAMELIST_ARRAYS:
            inp[key] = namelist_array(key, inp.get(key))
        inp["other_input"] = as_list(inp.get("other_input"))
        opt["inp"] = job_inp(inp, opt["model"]["data"])
    inp_content = template.render(opt)
    workdir = opt["model"]["data"]["workdir"]
    inp_script = run_file(workdir, run_name, segment_number, "inp")
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    with open(inp_script, "w", encoding="utf-8") as f:
        f.write(inp_content)    
    return inp_script

def _secondary_window(inp):
    return [parse_run_datetime(inp[key], key) for key in ("secondary_start_time", "secondary_stop_time")]


def coupled_handoff_regrids(engage_options):
    """True when the coupled grid differs from the spin-up grid, so the hand-off must regrid."""
    horires_standalone = float(engage_options["horires"])
    horires_coupled = float(engage_options["horires_coupled"])
    return ((horires_standalone, resolution_solver(horires_standalone)[0])
            != (horires_coupled, resolution_solver(horires_coupled)[0]))


def handoff_job_lines(machine, conda_env, script):
    """Job lines that activate conda_env (machines.yaml conda_activate) and run script; a batch
    job does not initialise conda itself."""
    activate = JobGen(CONFIG_DIR, machine=machine).m.get("conda_activate")
    if not activate:
        raise ValueError(f"config/machines.yaml has no conda_activate lines for machine {machine!r}: "
                         f"cannot write the job lines that run {script}.")
    return [line.format(conda_env=conda_env) for line in as_list(activate)] + [f"python {script}"]


def segment_plan(options, run_name, first_source=None):
    """(plans, pristart_times, pristop_times, last_prim_file) for a segmented run, writing
    nothing; every per-segment refusal is raised here."""
    segment = parse_segment(options["inp"].get("segment"))
    if segment is None:
        raise ValueError(f"A segmented run needs a segment length inp.segment = 'D H M S' "
                         f"(e.g. '5 0 0 0'); got {options['inp'].get('segment')!r}.")
    segment_times = segment_time(options["inp"]["start_time"], options["inp"]["stop_time"], segment)
    pri_files = 0
    sec_files = 0
    og_options = copy.deepcopy(options)
    PRIHIST = og_options["inp"]["PRIHIST"]
    MXHIST_PRIM = og_options["inp"]["MXHIST_PRIM"]
    SECHIST = og_options["inp"]["SECHIST"]
    MXHIST_SECH = og_options["inp"]["MXHIST_SECH"]
    PRIHIST_split = [int(i) for i in PRIHIST.split()]
    SECHIST_split = [int(i) for i in SECHIST.split()]
    step = _run_step(og_options)
    sec_window_start, sec_window_stop = _secondary_window(og_options["inp"])
    # histdir for the tool's own paths; histdir_job as the namelist names it.
    histdir = og_options["model"]["data"]["histdir"]
    histdir_job = history_dir(og_options["model"]["data"])
    job_name = og_options["simulation"]["job_name"]
    plans = []
    pristart_times = []
    pristop_times = []
    previous_START_YEAR = previous_PRISTOP = None
    for segment_number, segment in enumerate(segment_times):
        segment_options = copy.deepcopy(og_options)
        segment_start = segment[0]
        segment_stop = segment[1]
        segment_options["simulation"]["job_name"] =  job_name+"-{:02d}".format(segment_number +1)
        segment_START_YEAR, segment_START_DAY, segment_PRISTART, segment_PRISTOP = inp_pri_date(segment_start,segment_stop)
        if segment_number == 0:
            segment_options["inp"]["SOURCE"] = og_options["inp"]["SOURCE"] if first_source is None else first_source
            segment_options["inp"]["SOURCE_START"] = og_options["inp"]["SOURCE_START"]
        else:
            # Continue from the previous segment's last history, with the mtime nc_wrhist stored.
            segment_options["inp"]["SOURCE"] = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
            segment_options["inp"]["SOURCE_START"] = ' '.join(map(str, nc_wrhist_mtime(previous_PRISTOP, previous_START_YEAR)))
        segment_options["inp"]["START_YEAR"] = segment_START_YEAR
        segment_options["inp"]["START_DAY"] = segment_START_DAY
        segment_options["inp"]["PRISTART"] = ' '.join(map(str, segment_PRISTART))
        segment_options["inp"]["PRISTOP"] = ' '.join(map(str, segment_PRISTOP))
        previous_START_YEAR, previous_PRISTOP = segment_START_YEAR, segment_PRISTOP

        # PRIHIST and MXHIST_PRIM are fitted to each segment, which may be shorter than a full one.
        span_sec = dhms_to_seconds(segment_PRISTOP) - dhms_to_seconds(segment_PRISTART)
        segment_PRIHIST = fit_cadence(PRIHIST_split, span_sec, step)
        segment_MXHIST_PRIM = segment_mxhist(span_sec // dhms_to_seconds(segment_PRIHIST), MXHIST_PRIM)
        segment_OUTPUT, pri_files = inp_pri_out(segment_start, segment_stop, segment_PRIHIST, segment_MXHIST_PRIM, pri_files, histdir_job, run_name)
        segment_options["inp"]["PRIHIST"] = ' '.join(map(str, segment_PRIHIST))
        if segment_MXHIST_PRIM != int(MXHIST_PRIM):
            segment_options["inp"]["MXHIST_PRIM"] = segment_MXHIST_PRIM
        segment_options["inp"]["OUTPUT"] = segment_OUTPUT

        # Secondary history covers the part of the secondary window inside this segment; outside
        # it every secondary key is omitted (input.F requires all of them once any is given).
        sec_start = max(parse_run_datetime(segment_start, "segment start"), sec_window_start)
        sec_stop = min(parse_run_datetime(segment_stop, "segment stop"), sec_window_stop)
        sec_span_sec = int((sec_stop - sec_start).total_seconds())
        if sec_span_sec <= 0:
            for key in ("SECSTART", "SECSTOP", "SECHIST", "SECOUT", "MXHIST_SECH"):
                segment_options["inp"][key] = None
            segment_options["inp"]["SECFLDS"] = []
        else:
            sec_start_str = sec_start.strftime(RUN_DATETIME_FORMAT)
            sec_stop_str = sec_stop.strftime(RUN_DATETIME_FORMAT)
            segment_SECHIST = fit_cadence(SECHIST_split, sec_span_sec, step)
            segment_SECSTART, segment_SECSTOP = inp_sec_date(sec_start_str, sec_stop_str, segment_SECHIST)
            segment_options["inp"]["SECSTART"] = ' '.join(map(str, segment_SECSTART))
            segment_options["inp"]["SECSTOP"] = ' '.join(map(str, segment_SECSTOP))
            if segment_SECHIST != SECHIST_split:
                segment_options["inp"]["SECHIST"] = ' '.join(map(str, segment_SECHIST))
            segment_SECOUT, sec_files = inp_sec_out(sec_start_str, sec_stop_str, segment_SECHIST, MXHIST_SECH, sec_files, histdir_job, run_name)
            segment_options["inp"]["SECOUT"] = segment_SECOUT
        segment_options["model"]["data"]["log_file"] = os.path.join( options["model"]["data"]["workdir"],f"{run_name}-{'{:02d}'.format(segment_number+1)}.out")
        plans.append(segment_options)
        pristart_times.append(segment_options["inp"]["PRISTART"])
        pristop_times.append(segment_options["inp"]["PRISTOP"])
    last_prim_file = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files)}.nc"
    return plans, pristart_times, pristop_times, last_prim_file


def _run_step(options):
    """inp.STEP in seconds, else the resolution default."""
    step = options.get("inp", {}).get("STEP")
    if step is None or str(step).strip().lower() in ("", "none", "null"):
        return int(resolution_solver(options["model"]["specification"]["horires"])[3])
    return int(float(step))


def _handoff_script(segment_options, engage_options, last_prim_file):
    """(path, text) of the script that regrids the spin-up's last history onto the coupled grid."""
    data = segment_options["model"]["data"]
    path = os.path.join(data["workdir"], 'tiegcm_resolution_upscale.py')
    horires_coupled = engage_options["horires_coupled"]
    vertres_coupled = resolution_solver(horires_coupled, engage_options)[0]
    SOURCE_coupling = os.path.join(engage_options["parentdir"], f'{engage_options["job_name"]}_prim.nc')
    # Paths in the script are relative to its own directory.
    here = lambda p: job_path(p, data, dot=False)      # noqa: E731
    text = ("import os\n"
            "import sys\n"
            "here = os.path.dirname(os.path.abspath(__file__))\n"
            "sys.path.append(os.path.join(os.environ['TIEGCMHOME'], 'tiegcmrun'))\n"
            "import interpolation\n"
            f"interpolation.interpic(os.path.join(here, '{here(last_prim_file)}'),{float(horires_coupled)},{float(vertres_coupled)},{float(segment_options['model']['specification']['zitop'])},os.path.join(here, '{here(SOURCE_coupling)}'))\n")
    return path, text


def segment_jobs(options, run_name, pbs, engage_options=None, first_source=None):
    """Plan a segmented run and render its job scripts, writing nothing.

    Returns (plans, jobs, pristart_times, pristop_times, last_prim_file); jobs[i] is
    (pbs text or None, (hand-off path, text) or None)."""
    plans, pristart_times, pristop_times, last_prim_file = segment_plan(options, run_name, first_source)
    last_segment_time = len(plans) - 1
    jobs = []
    for segment_number, segment_options in enumerate(plans):
        data = segment_options["model"]["data"]
        data["input_file"] = run_file(data["workdir"], run_name, segment_number, "inp")
        pbs_text = handoff = None
        # On linux only an engage spin-up gets job scripts; a standalone run gets one run script.
        if pbs == True and (options["simulation"]["hpc_system"] != "linux" or engage_options is not None):
            if segment_number == last_segment_time and engage_options != None and coupled_handoff_regrids(engage_options):
                handoff = _handoff_script(segment_options, engage_options, last_prim_file)
                segment_options["job"]["job_chain"] = handoff_job_lines(
                    options["simulation"]["hpc_system"], engage_options["conda_env"],
                    job_path(handoff[0], data))
            pbs_text = render_pbs(segment_options)
        jobs.append((pbs_text, handoff))
    return plans, jobs, pristart_times, pristop_times, last_prim_file


def write_segment_jobs(planned, run_name):
    """Write what segment_jobs planned.

    Returns (inp_files, pbs_files, log_files, pristart_times, pristop_times, last_prim_file)."""
    plans, jobs, pristart_times, pristop_times, last_prim_file = planned
    inp_files = []
    pbs_files = []
    log_files = []
    for segment_number, segment_options in enumerate(plans):
        pbs_text, handoff = jobs[segment_number]
        segment_options["model"]["data"]["input_file"] = create_inp_scripts(segment_options,run_name,segment_number)
        pbs_script = None
        if handoff is not None:
            with open(handoff[0], "w", encoding="utf-8") as f:
                f.write(handoff[1])
        if pbs_text is not None:
            pbs_script = create_pbs_scripts(segment_options, run_name, segment_number, content=pbs_text)
        inp_files.append(segment_options["model"]["data"]["input_file"])
        pbs_files.append(pbs_script)
        log_files.append(segment_options["model"]["data"]["log_file"])
    return inp_files, pbs_files,log_files, pristart_times, pristop_times, last_prim_file


def segment_inp_pbs(options, run_name, pbs, engage_options=None, first_source=None):
    """Plan, then write, a segmented run's .inp (and .pbs) files; see write_segment_jobs.

    first_source overrides the first segment's SOURCE. Nothing is written if any segment fails."""
    return write_segment_jobs(segment_jobs(options, run_name, pbs, engage_options, first_source), run_name)
