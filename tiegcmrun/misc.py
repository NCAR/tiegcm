"""Utility functions for tiegcmrun: time segmentation, resolution and resource defaults, data-file
lookup, executable inspection and option-value checks."""

import os
import re
import copy
import json
import fnmatch
import argparse
from functools import lru_cache
from datetime import datetime, timedelta
import xarray as xr
from numpy import pad

from jobgen import spell

ENV_HINT = "source tiegcm/tiegcmrun/setEnvironment.sh (sets TIEGCMHOME and TIEGCMDATA)"


def tiegcm_env(name):
    """The value of an environment variable, or None when unset (read at call time so -h works)."""
    return os.environ.get(name) or None


def require_env(name, hint=ENV_HINT):
    """Return os.environ[name], or exit with a one-line message when it is unset."""
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"ERROR: environment variable {name} is not set; {hint}.")
    return value


def as_list(value, tokens=False):
    """Normalise a list-valued option to a list; [None] means no entries.

    Plain text is split on newlines and ';' (never on spaces); tokens=True also splits on commas
    and whitespace."""
    import ast
    if isinstance(value, (list, tuple)):
        return list(value)
    if value is None:
        return [None]
    s = str(value).strip()
    if s.lower() in ("", "none", "null", "[none]", "[null]"):
        return [None]
    if s.startswith("[") and s.endswith("]"):
        for parse in (ast.literal_eval, json.loads):
            try:
                parsed = parse(s)
            except (ValueError, SyntaxError, TypeError):
                continue
            if isinstance(parsed, (list, tuple)):
                return list(parsed) or [None]
    parts = []
    for line in s.replace(";", "\n").splitlines():
        line = line.strip()
        # unquote a single quoted entry, but not 'a' b 'c'
        if len(line) >= 2 and line[0] == line[-1] and line[0] in "'\"" and line[0] not in line[1:-1]:
            line = line[1:-1].strip()
        if tokens:
            parts += [w.strip("'\"") for w in line.replace(",", " ").split()]
        elif line:
            parts.append(line)
    return [p for p in parts if p] or [None]

def get_mtime(file_path):
    """The file's mtime records, each padded to 4 entries."""
    with xr.open_dataset(file_path) as ds:
        if 'mtime' not in ds.variables:
            raise ValueError(f"{file_path} has no 'mtime' variable")
        mtime_data = ds['mtime'].values
    mtime_arr = pad(mtime_data, [(0, 0), (0, max(4 - mtime_data.shape[1], 0))], mode='constant').tolist()
    return mtime_arr

def segment_time(start_time_str, stop_time_str, interval_array):
    """Split [start, stop] into [[start, end], ...] steps of interval_array ([d, h, m, s]),
    breaking every segment at Jan 1."""
    start = datetime.strptime(start_time_str, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time_str, '%Y-%m-%dT%H:%M:%S')

    days, hours, minutes, seconds = interval_array

    delta = timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)

    intervals = []
    current = start
    temp_delta = 0
    while current < stop:
        if temp_delta != 0:
            next_time = min(current + temp_delta, stop)
            temp_delta = 0
        else:
            next_time = min(current + delta, stop)
        if current.year != next_time.year:
            if next_time.month != 1 or next_time.day != 1 or next_time.hour != 0 or next_time.minute != 0 or next_time.second != 0:
                year_boundary = datetime(current.year + 1, 1, 1, 0, 0, 0)
                intervals.append([
                    current.strftime('%Y-%m-%dT%H:%M:%S'),
                    year_boundary.strftime('%Y-%m-%dT%H:%M:%S')
                ])
                current = year_boundary
                temp_delta = next_time - year_boundary
            else:
                intervals.append([
                current.strftime('%Y-%m-%dT%H:%M:%S'),
                next_time.strftime('%Y-%m-%dT%H:%M:%S')
                ])
                current = next_time
        else:
            intervals.append([
                current.strftime('%Y-%m-%dT%H:%M:%S'),
                next_time.strftime('%Y-%m-%dT%H:%M:%S')
            ])
            current = next_time
    return intervals

def valid_bench(value):
    """argparse type: a benchmark name from benchmarks.json."""
    if value is None:
        return value
    benchmarks_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "benchmarks.json")
    with open(benchmarks_file, "r", encoding="utf-8") as f:
        valid = [k for k in json.load(f) if k != "version"]
    if value not in valid:
        raise argparse.ArgumentTypeError(
            f"{value} is not a valid benchmark option. Valid options: {', '.join(valid)}")
    return value

# Magnetic grid resolution (deg) -> defs.h NRES_GRID.
MRES_TO_NRES_GRID = {2.0: 5, 1.0: 6, 0.5: 7}

MAKE_FRAGMENTS = {"derecho": "Make.intel_de", "aitken": "Make.intel_at", "linux": "Make.intel_linux"}


def mres_to_nres_grid(mres):
    try:
        return MRES_TO_NRES_GRID[float(mres)]
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"Unsupported magnetic grid resolution: {mres!r} "
                         f"(expected one of 2, 1, 0.5)") from None


def in_execdir(execdir, exe):
    """exe joined to execdir, unless it already lies under execdir."""
    if os.path.normpath(exe).startswith(os.path.normpath(execdir) + os.sep):
        return exe
    return os.path.join(execdir, exe)


def default_make_file(modeldir, hpc_system):
    """<modeldir>/scripts/Make.* for hpc_system, or None."""
    frag = MAKE_FRAGMENTS.get(str(hpc_system).strip().lower())
    return os.path.join(modeldir, "scripts", frag) if frag else None


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
    else:
        raise ValueError(f"Unsupported horizontal resolution: {horires!r}")

    nres_grid = mres_to_nres_grid(mres)
    
    if engage_options != None:
        STEP = engage_options["STEP"]
    
    return vertres, mres, nres_grid, STEP

# Seasonal SOURCE file prefixes and the day of year of the history each holds.
SEASONAL_SOURCES = (("mareqx", 81), ("junsol", 173), ("sepeqx", 265), ("decsol", 356))


def nearest_season(day_of_year):
    """The seasonal SOURCE prefix nearest to day_of_year, wrapping the year; ties go to the
    earlier season."""
    best = None
    for name, day in SEASONAL_SOURCES:
        forward = (int(day_of_year) - day) % 365       # days since this season's history
        key = (min(forward, 365 - forward), forward)
        if best is None or key < best[0]:
            best = (key, name)
    return best[1]


def select_source_defaults(options, option_descriptions, data_dir=None):
    """The default inp.SOURCE: the nearest seasonal file at the flux level's F10.7, or None."""
    start_time = options["inp"]["start_time"]
    time_dhms = time_to_dhms(start_time)
    flux_level = options["inp"]["solar_flux_level"]
    if flux_level == "low":
        f107 = 70
    elif flux_level == "medium":
        f107 = 140
    elif flux_level == "high":
        f107 = 200
    else:
        raise ValueError(f"Unknown solar_flux_level: {flux_level!r} (expected low|medium|high)")
    return find_file(f"{nearest_season(time_dhms[0])}_f{f107}*", data_dir or tiegcm_env("TIEGCMDATA"))

# TIE-GCM task-count rules, ported from mpi.F and input.F mkntask: ntask must be 1, 2 or not prime,
# factor onto both grids, and leave every task at least 4 longitudes.

def tiegcm_grid(horires, nres_grid=None):
    """(nlonp4, nlat, nmlonp1, nmlat) of the model grid (params.F)."""
    horires = float(horires)
    if nres_grid in (None, "", "None"):
        nres_grid = resolution_solver(horires)[2]
    nres_grid = int(nres_grid)
    nlat, nlon = int(180 / horires), int(360 / horires)
    return nlon + 4, nlat, 5 * 2 ** (nres_grid - 1) + 1, 3 * 2 ** nres_grid + 1


def mkntask(ntask, nx, ny):
    """input.F mkntask: the (ntaskx, ntasky) the model picks, or None when none fits."""
    choices = [(i, j) for i in range(1, nx + 1) if ntask % i == 0
               for j in (ntask // i,) if j <= ny]
    if not choices:
        return None
    square = [c for c in choices if c[0] == c[1]]
    if square:
        return square[-1]
    best, gap = None, nx * ny
    for i, j in choices:
        if abs(i - j) < gap:
            best, gap = (i, j), abs(i - j)
    x, y = best
    return (y, x) if x < y else (x, y)


def ntask_problem(ntask, horires, nres_grid=None):
    """Why TIE-GCM refuses ntask MPI tasks, or None when it accepts them."""
    try:
        ntask = int(ntask)
    except (TypeError, ValueError):
        return f"{ntask!r} is not a whole number of tasks"
    if ntask < 1:
        return "fewer than 1 task"
    if ntask not in (1, 2) and all(ntask % i for i in range(2, int(ntask ** 0.5) + 1)):
        return "a prime task count"
    nlonp4, nlat, nmlonp1, nmlat = tiegcm_grid(horires, nres_grid)
    for label, nx, ny in (("geographic", nlonp4, nlat), ("magnetic", nmlonp1, nmlat)):
        pick = mkntask(ntask, nx, ny)
        if pick is None:
            return f"no {label} task layout i x j = {ntask} with i <= {nx}, j <= {ny}"
        if nx // pick[0] < 4:
            return (f"the {nx} {label} longitudes split over {pick[0]} tasks give fewer than 4 "
                    f"per task")
    return None


def nearest_valid_ntasks(ntask, horires, nres_grid=None):
    """(nearest valid task count below, nearest above); None where none exists."""
    ntask = int(ntask)
    below = next((n for n in range(ntask - 1, 0, -1) if ntask_problem(n, horires, nres_grid) is None), None)
    nlonp4, nlat, _, _ = tiegcm_grid(horires, nres_grid)
    above = next((n for n in range(ntask + 1, nlonp4 * nlat + 1)
                  if ntask_problem(n, horires, nres_grid) is None), None)
    return below, above


def ntask_error(ntask, horires, nres_grid=None, where="job.nprocs"):
    """An error message for an invalid task count, or None when it is valid."""
    problem = ntask_problem(ntask, horires, nres_grid)
    if problem is None:
        return None
    try:
        below, above = nearest_valid_ntasks(ntask, horires, nres_grid)
        near = " or ".join(str(n) for n in (below, above) if n is not None) or "none"
    except (TypeError, ValueError):
        near = "none"
    return (f"{where} = {ntask} MPI tasks: TIE-GCM at horires {float(horires)} refuses it ({problem}). "
            f"Nearest task counts it accepts: {near} (select * mpiprocs).")


# The aitken (NAS) node models: :model= tag -> (cores per node, default MPI ranks per node).
AITKEN_NODE_MODELS = {"mil_ait": (128, 128), "rom_ait": (128, 128), "cas_ait": (40, 36),
                      "bro": (28, 24), "has": (24, 24), "ivy": (20, 18), "san": (16, 12),
                      "sky_ele": (40, 36), "bro_ele": (28, 24)}


def _feasible_layout(select, mpiprocs, horires, nres_grid):
    """The largest (select, mpiprocs) not above the given one whose product TIE-GCM accepts."""
    for nodes in range(int(select), 0, -1):
        for ranks in range(int(mpiprocs), 0, -1):
            if ntask_problem(nodes * ranks, horires, nres_grid) is None:
                return nodes, ranks
    return 1, 1


def select_resource_defaults(options, option_descriptions, queue_sized=True):
    """Default (select, ncpus, mpiprocs) for the machine, capped by the queue and reduced to a
    task count TIE-GCM accepts."""
    spec = options["model"]["specification"]
    horires = spec["horires"]
    nres_grid = spec.get("nres_grid")
    hpc_platform = options["simulation"]["hpc_system"]
    if hpc_platform == "derecho":
        select_default, ncpus_default, mpiprocs_default = 3, 128, 96
    elif hpc_platform == "aitken":
        model = ((options.get("job") or {}).get("resource") or {}).get("model")
        if model not in AITKEN_NODE_MODELS:
            raise ValueError(f"job.resource.model {model!r} is not an aitken node model "
                             f"({', '.join(AITKEN_NODE_MODELS)}).")
        ncpus_default, mpiprocs_default = AITKEN_NODE_MODELS[model]
        select_default = 1
    else:
        from jobgen import JobGen
        cores = int(JobGen(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config"),
                           machine=hpc_platform).m["cores_per_node"])
        select_default, ncpus_default, mpiprocs_default = 3, cores, min(96, cores)
    queue = (options.get("job") or {}).get("queue")
    if queue_sized and not is_unset(queue):
        from output_solver import queue_limits
        ncpus_max = queue_limits(hpc_platform, str(queue).strip()).get("ncpus_max")
        if ncpus_max is not None:
            select_default = min(select_default, max(1, int(ncpus_max) // int(ncpus_default)))
    select_default, mpiprocs_default = _feasible_layout(select_default, mpiprocs_default,
                                                        horires, nres_grid)
    return select_default, ncpus_default, mpiprocs_default

@lru_cache(maxsize=None)
def _walk_files(path):
    """Cached [(name, fullpath), ...] listing of path, in os.walk order."""
    listing = []
    for root, dirs, files in os.walk(path):
        for name in files:
            listing.append((name, os.path.join(root, name)))
    return listing


GSWM_PATTERNS = {
    "GSWM_MI_DI_NCFILE":  "*gswm_diurn_{horires}d_99km*",
    "GSWM_MI_SDI_NCFILE": "*gswm_semi_{horires}d_99km*",
    "GSWM_NM_DI_NCFILE":  "*gswm_nonmig_diurn_{horires}d_99km*",
    "GSWM_NM_SDI_NCFILE": "*gswm_nonmig_semi_{horires}d_99km*",
}
# Non-migrating GSWM tides: input.F shuts down on them at 5 degrees.
GSWM_NONMIGRATING = ("GSWM_NM_DI_NCFILE", "GSWM_NM_SDI_NCFILE")


def gswm_allowed(key, horires):
    return not (key in GSWM_NONMIGRATING and float(horires) == 5.0)


def he_coefs_pattern(horires):
    """The helium-coefficient file input.F defaults to on this grid."""
    return "*he_coefs_sres*" if float(horires) == 5.0 else "*he_coefs_dres*"


# &tgcm_input array keys (input.F): name -> (max entries, element type). Numeric arrays are
# rendered unquoted (a quoted value fails the read with iostat 64).
NAMELIST_NUMERIC_ARRAYS = {
    "TIDE": (10, float), "TIDE2": (2, float), "NUDGE_REFDATE": (2, int),
    "NUDGE_SPONGE": (2, float), "NUDGE_DELTA": (2, float), "NUDGE_POWER": (2, float),
}
NAMELIST_STRING_ARRAYS = {"NUDGE_FLDS": (500, str), "NUDGE_NCFILE": (500, str)}
NAMELIST_BOOL_ARRAYS = {"NUDGE_LBC": 500, "NUDGE_F4D": 500}
NAMELIST_ARRAYS = {**NAMELIST_NUMERIC_ARRAYS, **NAMELIST_STRING_ARRAYS,
                   **{k: (n, bool) for k, n in NAMELIST_BOOL_ARRAYS.items()}}


def namelist_array(name, value):
    """A namelist array value as the list template.inp renders, or None when unset.

    Numeric entries are validated but kept as the user's text so a replay renders the same bytes;
    logical entries become .true./.false.; raises ValueError naming the key."""
    size, kind = NAMELIST_ARRAYS[name]
    items = [x for x in as_list(value, tokens=True) if x is not None]
    items = [t for x in items for t in re.split(r"[,\s]+", str(x).strip())]
    items = [x.strip("'\"") for x in items]
    items = [x for x in items if x != ""]
    if not items:
        return None
    if len(items) > size:
        raise ValueError(f"{name} takes at most {size} values; got {len(items)}")
    if kind is bool:
        out = []
        for x in items:
            try:
                out.append(".true." if _boolean(name, x) else ".false.")
            except ValueError:
                raise ValueError(f"{name} entry {x!r} is not true or false") from None
        return out
    for x in items:
        try:
            if kind is int:
                int(x)
            elif kind is float:
                float(x.lower().replace("d", "e"))
        except ValueError:
            raise ValueError(f"{name} entry {x!r} is not {'an integer' if kind is int else 'a number'}")
    return items


# Fallback &tgcm_input member list, used only when input.F cannot be read.
_TGCM_INPUT_MEMBERS = (
    "label step nstep_sub source source_start output pristart pristop prihist secout secstart "
    "secstop sechist secflds potential_model eddy_dif dynamo tide tide2 tide3m3 f107 f107a power "
    "ctpoten bximf byimf bzimf swvel swden al kp colfac joulefac aurora gpi_ncfile "
    "gswm_mi_di_ncfile gswm_mi_sdi_ncfile gswm_nm_di_ncfile gswm_nm_sdi_ncfile mxhist_prim "
    "mxhist_sech ntask_lat ntask_lon start_day start_year calendar_advance see_ncfile "
    "ctpoten_time power_time bximf_time byimf_time bzimf_time kp_time al_time swden_time "
    "swvel_time indices_interp imf_ncfile saber_ncfile tidi_ncfile sech_nbyte f107_time "
    "f107a_time hpss_path current_pg current_kq calc_helium bgrddata_ncfile ctmt_ncfile "
    "electron_heating duff amienh amiesh amie_ibkg he_coefs_ncfile enforce_n2 et saps "
    "subaur_data doeclipse eclipse_list oneway mixfile nudge_ncpre nudge_ncfile nudge_ncpost "
    "nudge_flds nudge_lbc nudge_f4d nudge_use_refdate nudge_refdate nudge_sponge nudge_delta "
    "nudge_power nudge_alpha opdiffcap opdiffrate opdifflev opfloor oprate oplev oplatwidth "
    "te_cap ti_cap").split()


def _parse_tgcm_input(path):
    """Lower-case member names of namelist/tgcm_input/ in input.F, or None."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        m = re.match(r"^\s+namelist\s*/\s*tgcm_input\s*/(.*)$", line, re.I)
        if not m:
            continue
        body = [m.group(1)]
        j = i + 1
        while j < len(lines) and re.match(r"^     [^ 0]", lines[j]):   # continuation lines
            body.append(lines[j][6:])
            j += 1
        text = " ".join(b.split("!")[0] for b in body)
        names = {n.strip().lower() for n in text.split(",") if n.strip()}
        return names or None
    return None


def namelist_members(modeldir=None):
    """&tgcm_input member names from <modeldir> or $TIEGCMHOME src/input.F, else the fallback."""
    for root in (modeldir, tiegcm_env("TIEGCMHOME")):
        if root:
            names = _parse_tgcm_input(os.path.join(os.path.expanduser(root), "src", "input.F"))
            if names:
                return names
    return set(_TGCM_INPUT_MEMBERS)


def _unquoted(line):
    return re.sub(r"'[^']*'|\"[^\"]*\"", "''", line)


def other_input_keys(line):
    """The namelist keys an other_input line assigns."""
    s = _unquoted(str(line)).split("!")[0]
    return [m.group(1) for m in re.finditer(r"(?:^|[,\s])([A-Za-z_]\w*)\s*(?:\([^)]*\))?\s*=", s)]


def unknown_other_input_keys(lines, modeldir=None):
    """Upper-case other_input keys that are not &tgcm_input members (the namelist read fails)."""
    members = namelist_members(modeldir)
    bad = []
    for line in lines or []:
        if line is None:
            continue
        for key in other_input_keys(line):
            if key.lower() not in members and key.upper() not in bad:
                bad.append(key.upper())
    return bad


def other_input_error(keys):
    return (f"other_input: {', '.join(keys)} {'is not a' if len(keys) == 1 else 'are not'} "
            f"&tgcm_input namelist member{'s' if len(keys) > 1 else ''}: the model would stop "
            f"reading {'it' if len(keys) == 1 else 'them'}. Check the spelling.")


# SECFLDS names are character(len=16), at most 500 (params.F); the model compares them
# case-sensitively and stops on a duplicate.
SECFLDS_NAME_LEN = 16
SECFLDS_MAX = 500
SECFLDS_MANDATORY = ("TN", "O2", "O1", "Z", "ZG", "ZMAG")
_CURLY_QUOTES = str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"'})
# Fallback field list (gen_secflds_fields.py), used only when no model source can be read.
SECFLDS_FIELDS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "secflds_fields.json")


def ascii_quotes(text):
    return str(text).translate(_CURLY_QUOTES)


def scan_secflds_fields(src):
    """Every secondary-history field name found in the model source tree src, or None."""
    import glob
    files = sorted(glob.glob(os.path.join(src, "*.F")) + glob.glob(os.path.join(src, "*.F90")))
    if not files:
        return None
    names = set(SECFLDS_MANDATORY)
    patterns = (re.compile(r"short_name\s*=\s*\(?\s*[\"']([A-Za-z0-9_]+)[\"']"),
                re.compile(r"call\s+addfld\s*\(\s*'([A-Za-z0-9_]+)'", re.I),
                re.compile(r"call\s+mkdiag_\w+\s*\(\s*'([A-Za-z0-9_]+)'", re.I))
    # minor.F builds names from the species passed by comp_*.F: every template x every species.
    templates = re.compile(r"call\s+addfld\s*\(\s*(?:'([A-Za-z0-9_]*)'\s*//\s*)?name\s*"
                           r"(?://\s*'([A-Za-z0-9_]*)')?\s*,", re.I)
    minor_call = re.compile(r"call\s+minor\s*\([^']*?'([A-Za-z0-9_]+)'\s*\)", re.I | re.S)
    minor_templates, minor_species = set(), set()
    for path in files:
        try:
            with open(path, errors="replace") as f:
                lines = f.readlines()
        except OSError:
            continue
        active = []
        for line in lines:
            s = line.lstrip()
            if s.startswith("!") or (line[:1] in "cC*" and path.endswith(".F")):
                continue          # fixed-form comment: c/C/* in column 1
            active.append(s)
            for pattern in patterns:
                names.update(m.group(1) for m in pattern.finditer(s))
        if os.path.basename(path).lower() == "minor.f":
            for s in active:
                minor_templates.update((m.group(1) or "", m.group(2) or "") for m in templates.finditer(s))
        text = "".join(re.sub(r"!.*", "", s) for s in active)
        minor_species.update(m.group(1) for m in minor_call.finditer(text))
    names.update(f"{pre}{sp}{post}" for pre, post in minor_templates for sp in minor_species)
    return names


@lru_cache(maxsize=None)
def _secflds_fields_cached(src):
    names = scan_secflds_fields(src) if src else None
    return frozenset(names) if names else None


def secflds_fields(modeldir=None):
    """Secondary-history field names (exact case) from <modeldir> or $TIEGCMHOME src, else the
    shipped secflds_fields.json."""
    for root in (modeldir, tiegcm_env("TIEGCMHOME")):
        if root:
            names = _secflds_fields_cached(os.path.join(os.path.expanduser(str(root)), "src"))
            if names:
                return names
    try:
        with open(SECFLDS_FIELDS_FILE, encoding="utf-8") as f:
            return frozenset(json.load(f)["fields"])
    except (OSError, ValueError, KeyError):
        return frozenset(SECFLDS_MANDATORY)


def secflds_list(value, modeldir=None, where="SECFLDS", fix="remove or correct them", short=False):
    """(field names for the deck, normalisation notes) for a SECFLDS value.

    Fixes case where unambiguous and drops repeats; raises ValueError for names that are too long,
    unknown, or too many."""
    fields = secflds_fields(modeldir)
    by_lower = {}
    for name in fields:
        by_lower.setdefault(name.lower(), []).append(name)
    items = value if isinstance(value, (list, tuple)) else as_list(ascii_quotes(value), tokens=True)
    names, notes, too_long, unknown, ambiguous = [], [], [], [], []
    for raw in items:
        if raw is None:
            continue
        text = ascii_quotes(raw)
        words = [w.strip("'\"") for w in text.replace(",", " ").split()]
        for word in (w for w in words if w):
            if word != str(raw).strip():
                notes.append(f"{str(raw).strip()!r} -> {word}")
            if len(word) > SECFLDS_NAME_LEN:
                too_long.append(word)
                continue
            if word not in fields:
                matches = by_lower.get(word.lower(), [])
                if len(matches) == 1:
                    notes.append(f"{word} -> {matches[0]}")
                    word = matches[0]
                elif matches:
                    ambiguous.append(f"{word} ({' or '.join(sorted(matches))})")
                    continue
                else:
                    unknown.append(word)
                    continue
            if word in names:
                notes.append(f"repeated {word} dropped")
                continue
            names.append(word)
    errors = []
    if too_long:
        errors.append(f"{', '.join(too_long)}: longer than {SECFLDS_NAME_LEN} characters (TIE-GCM "
                      f"cuts field names to {SECFLDS_NAME_LEN})")
    if ambiguous:
        errors.append(f"{', '.join(ambiguous)}: the case matters (two different fields)")
    if unknown:
        import difflib
        hints = []
        for word in unknown:
            close = difflib.get_close_matches(word.upper(), sorted(fields, key=str.upper), n=2, cutoff=0.75)
            if close:
                hints.append(f"{word} -> {' or '.join(close)}")
        errors.append(f"{', '.join(unknown)}: {'is not a' if len(unknown) == 1 else 'are not'} "
                      f"TIE-GCM secondary-history field{'s' if len(unknown) > 1 else ''}"
                      + (f" (did you mean: {'; '.join(hints)})" if hints else ""))
    if len(names) > SECFLDS_MAX:
        errors.append(f"{len(names)} names: TIE-GCM takes at most {SECFLDS_MAX}")
    if errors and short:
        raise ValueError(f"{where}: " + "; ".join(errors) + "; the valid names are in "
                         "tiegcmrun/secflds_fields.json.")
    if errors:
        valid = ", ".join(sorted(fields, key=lambda n: (n.upper(), n)))
        bad = too_long + [a.split(" ")[0] for a in ambiguous] + unknown
        raise ValueError(f"{where}: " + "; ".join(errors) + f". To fix: {fix} ({', '.join(bad)})."
                         f" Valid SECFLDS names: {valid}")
    return (names or [None]), notes


def find_file(pattern, path):
    """The first file under path whose name matches pattern, or None."""
    for name, fullpath in _walk_files(path):
        if fnmatch.fnmatch(name, pattern):
            return fullpath
    return None

_YEARDAY_BOUNDS_RE = re.compile(r'_(\d+)-(\d+)\.nc$')

def _file_yearday_bounds(path):
    """(beg, end) YYYYDDD from a <prefix>_<beg>-<end>.nc filename, or None."""
    m = _YEARDAY_BOUNDS_RE.search(os.path.basename(path))
    return (int(m.group(1)), int(m.group(2))) if m else None

def select_latest_gpi(path):
    """The gpi_<beg>-<end>.nc file under path with the latest end date, else any gpi_* file."""
    latest, latest_end = None, -1
    for name, fullpath in _walk_files(path):
        if fnmatch.fnmatch(name, 'gpi_*.nc'):
            bounds = _file_yearday_bounds(fullpath)
            if bounds is not None and bounds[1] > latest_end:
                latest_end, latest = bounds[1], fullpath
    return latest if latest is not None else find_file('gpi_*', path)

def file_covers(file_path, start_time, stop_time):
    """True if the file's dated name covers the whole window [start_time, stop_time]."""
    if not file_path:
        return False
    bounds = _file_yearday_bounds(file_path)
    if bounds is None:
        return False
    beg, end = bounds
    def yd(t):
        d = datetime.strptime(t, '%Y-%m-%dT%H:%M:%S')
        return d.year * 1000 + d.timetuple().tm_yday
    return beg <= yd(start_time) and end >= yd(stop_time)

def benchmark_window(bench_inp):
    """(start_time, stop_time) ISO strings of a benchmarks.json 'inp' block."""
    year = int(bench_inp["START_YEAR"])
    def iso(mtime):
        d, h, m, s = (list(mtime) + [0, 0, 0, 0])[:4]
        t = datetime(year, 1, 1) + timedelta(days=int(d) - 1, hours=int(h), minutes=int(m),
                                             seconds=int(s))
        return t.strftime('%Y-%m-%dT%H:%M:%S')
    return iso(bench_inp["PRISTART"]), iso(bench_inp["PRISTOP"])

def resolve_benchmark_file(field, benchmark, bench_inp, path):
    """A benchmark's SOURCE / GPI_NCFILE / IMF_NCFILE file under path.

    Raises FileNotFoundError when it is missing or does not span the benchmark window."""
    raw = bench_inp.get(field)
    if field == "SOURCE":
        want = f"{benchmark}.nc"
        found = find_file(want, path)
    elif field == "GPI_NCFILE" and any(c in str(raw) for c in "*?["):
        want = raw
        found = select_latest_gpi(path)
    else:
        want = raw
        found = find_file(raw, path)
    if found is None:
        raise FileNotFoundError(
            f"benchmark {benchmark}: {field} file {want!r} was not found under TIEGCMDATA={path}. "
            f"Point TIEGCMDATA at the TIE-GCM 3.0 data tree that ships the benchmark files.")
    if field in ("GPI_NCFILE", "IMF_NCFILE") and _file_yearday_bounds(found) is not None:
        start, stop = benchmark_window(bench_inp)
        if not file_covers(found, start, stop):
            raise FileNotFoundError(
                f"benchmark {benchmark}: {field} file {found} does not span the benchmark window "
                f"{start}..{stop}.")
    return found

def _import_gcmprocpy(kind):
    try:
        import gcmprocpy
        return gcmprocpy
    except ImportError as e:
        raise RuntimeError(
            f"{kind} generation needs the 'gcmprocpy' package (>=1.5) in the active "
            f"environment. Install it (pip install gcmprocpy) or run tiegcmrun from a "
            f"conda env that has it."
        ) from e

def generate_gpi_file(start_time, stop_time, output_dir, verbose=True):
    """Write a GPI file for the window into output_dir and return its path.

    The 27-day F107a average is trailing, so it needs no data after stop_time."""
    gcmprocpy = _import_gcmprocpy("GPI")
    if not output_dir:
        raise RuntimeError("Cannot write a generated GPI file: the run work directory is not set.")
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    print(f"Generating GPI (27-day trailing F107a) for {start.date()} -> {stop.date()} via gcmprocpy ...")
    ds = gcmprocpy.generate_gpi(start=start, end=stop, window=27, centered=False, verbose=verbose)
    path = gcmprocpy.save_gpi(ds, output_dir=output_dir)
    print(f"Generated GPI file: {path}")
    return path

def generate_imf_file(start_time, stop_time, output_dir, verbose=True):
    """Write an OMNI solar-wind IMF file for the window into output_dir and return its path."""
    gcmprocpy = _import_gcmprocpy("IMF")
    if not output_dir:
        raise RuntimeError("Cannot write a generated IMF file: the run work directory is not set.")
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    print(f"Generating IMF (OMNI solar wind) for {start.date()} -> {stop.date()} via gcmprocpy ...")
    ds = gcmprocpy.generate_imf(start=start, end=stop, verbose=verbose)
    path = gcmprocpy.save_imf(ds, output_dir=output_dir)
    print(f"Generated IMF file: {path}")
    return path

def generate_pending(inp, workdir):
    """Replace any 'gen' GPI_NCFILE / IMF_NCFILE in inp with a generated file; returns inp."""
    start, stop = inp.get("start_time"), inp.get("stop_time")
    for key, gen in (("GPI_NCFILE", generate_gpi_file), ("IMF_NCFILE", generate_imf_file)):
        val = inp.get(key)
        if not (isinstance(val, str) and val.strip().lower() == "gen"):
            continue
        print(f"==> {key}: generating a file for {start}..{stop} into {workdir} ...")
        try:
            inp[key] = gen(start, stop, workdir)
        except Exception as e:
            raise RuntimeError(f"{key}: generation failed: {e}. Provide a file path instead, or run "
                               f"from an environment that has gcmprocpy (pip install gcmprocpy).") from e
        print(f"==> {key}: generated {inp[key]}")
        if not file_covers(inp[key], start, stop):
            print(f"WARNING {key}: generated file {inp[key]} does not span the full run "
                  f"({start}..{stop}); the indices source may lag real time.")
    return inp


def time_to_dhms(time_str):

    time = datetime.strptime(time_str, "%Y-%m-%dT%H:%M:%S")

    day = time.timetuple().tm_yday
    hour = time.hour
    minute = time.minute
    second = time.second
    
    return [day, hour, minute, second]

def seconds_to_dhms(seconds):
    days = seconds // (24 * 3600)
    seconds %= (24 * 3600)

    hours = seconds // 3600
    seconds %= 3600

    minutes = seconds // 60

    seconds %= 60
    
    return [days, hours, minutes, seconds]

def dhms_to_seconds(dhms):
    return dhms[0]*86400 + dhms[1]*3600 + dhms[2]*60 + dhms[3]

# An executable's grid is read from the nm -S sizes of the params_module arrays (8-byte reals);
# a coupled (-DGAMERA) build has mage_oneway_mp_* symbols.
ELF_MAGIC = b"\x7fELF"
TIEGCM_ZIBOT = -7.0
TIEGCM_ZITOP = 7.0
_NM_ARRAYS = {"nlon": "params_module_mp_glon_", "nlat": "params_module_mp_glat_",
              "nlevp1": "params_module_mp_zpmid_", "nmlat": "params_module_mp_gmlat_"}
_IDENTITY_CACHE = {}


def is_elf(path):
    try:
        with open(path, "rb") as f:
            return f.read(4) == ELF_MAGIC
    except (OSError, TypeError, ValueError):
        return False


def nm_signature(path):
    """{nlon, nlat, nlevp1, nmlat, mage_oneway_syms} from nm -S, or None when the symbol table
    cannot tell (nm fails, or the executable is stripped)."""
    import subprocess
    try:
        p = subprocess.run(["nm", "-S", str(path)], capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.SubprocessError):
        return None
    if p.returncode != 0:
        return None
    sizes, mage = {}, 0
    for line in p.stdout.splitlines():
        if "mage_oneway_mp_" in line:
            mage += 1
        f = line.split()
        if len(f) == 4:
            try:
                sizes[f[3]] = int(f[1], 16)
            except ValueError:
                pass
    out = {k: (sizes[s] // 8 if sizes.get(s) else None) for k, s in _NM_ARRAYS.items()}
    if not mage and all(v is None for v in out.values()):
        return None
    out["mage_oneway_syms"] = mage
    return out


def grid_from_nm(nlon, nlevp1, nmlat, zitop=None):
    """(horires, vertres, nres_grid) from the nm array sizes; ZIBOT is -7."""
    ztop = TIEGCM_ZITOP if zitop in (None, "", "None") else float(zitop)
    horires = 360 / nlon if nlon else None
    vertres = (ztop - TIEGCM_ZIBOT) / (nlevp1 - 1) if nlevp1 and nlevp1 > 1 else None
    nres_grid = None
    if nmlat and nmlat > 1 and (nmlat - 1) % 3 == 0:
        n = (nmlat - 1) // 3
        if n & (n - 1) == 0:                        # a power of two: nmlat = 3*2**nres_grid + 1
            nres_grid = n.bit_length() - 1
    return horires, vertres, nres_grid


def exe_identity(path, zitop=None):
    """A dict describing the executable at path (ELF, coupling, grid), read from the file; values
    that cannot be read are None. Never raises."""
    path = None if path in (None, "") else str(path)
    out = {"path": path, "exists": bool(path) and os.path.isfile(path), "elf": False,
           "size": None, "mage_oneway_syms": None, "coupling": None,
           "nlon": None, "nlat": None, "nlevp1": None, "nmlat": None,
           "horires": None, "vertres": None, "nres_grid": None}
    if not out["exists"]:
        return out
    try:
        st = os.stat(path)
    except OSError:
        out["exists"] = False
        return out
    key = (os.path.realpath(path), st.st_size, st.st_mtime_ns, zitop)
    if key in _IDENTITY_CACHE:
        cached = dict(_IDENTITY_CACHE[key])
        cached["path"] = path
        return cached
    out["size"] = st.st_size
    out["elf"] = is_elf(path)
    if out["elf"]:
        nm = nm_signature(path)
        if nm is not None:
            out.update({k: nm[k] for k in ("nlon", "nlat", "nlevp1", "nmlat", "mage_oneway_syms")})
            out["coupling"] = nm["mage_oneway_syms"] > 0
            out["horires"], out["vertres"], out["nres_grid"] = grid_from_nm(
                nm["nlon"], nm["nlevp1"], nm["nmlat"], zitop)
    _IDENTITY_CACHE[key] = dict(out)
    return out


def _first_line(path, limit=40):
    try:
        with open(path, "rb") as f:
            line = f.readline(limit).split(b"\n")[0]
        return line.decode("utf-8", "replace")
    except OSError:
        return ""


def exe_problems(identity, role, coupled, horires, vertres, nres_grid, zitop=None):
    """(not_elf message or None, [problems]) for an executable in this run: wrong build kind or
    a grid other than the run's. Nothing is reported for what cannot be read."""
    path = identity.get("path")
    if not identity.get("exists"):
        return None, []
    if not identity.get("elf"):
        return (f"{role} {path} is not an ELF executable ({identity.get('size')} bytes, first line "
                f"{_first_line(path)!r}): the job could not launch it as TIE-GCM"), []
    problems = []
    syms = identity.get("mage_oneway_syms")
    if syms is not None:
        if coupled and not syms:
            problems.append(f"{role} {path} is a standalone build (no mage_oneway_mp_* symbols; built "
                            f"without --coupling): as the coupled TIE-GCM it never exchanges with "
                            f"voltron, which then waits until the job's walltime")
        elif not coupled and syms:
            problems.append(f"{role} {path} is a coupled build ({syms} mage_oneway_mp_* symbols; "
                            f"built with --coupling): as a standalone TIE-GCM it waits for a voltron "
                            f"that is not there")
    def want(fn):
        try:
            return fn()
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    ztop = want(lambda: TIEGCM_ZITOP if zitop in (None, "", "None") else float(zitop))
    want_nlon = want(lambda: int(round(360 / float(horires))))
    want_nlevp1 = want(lambda: int(round((ztop - TIEGCM_ZIBOT) / float(vertres))) + 1)
    want_nmlat = want(lambda: 3 * 2 ** int(float(nres_grid)) + 1)
    got = []
    if want_nlon and identity.get("nlon") and identity["nlon"] != want_nlon:
        got.append(f"horires {identity['horires']:g} deg (nlon {identity['nlon']}), the run needs "
                   f"{float(horires):g} deg (nlon {want_nlon})")
    if want_nlevp1 and identity.get("nlevp1") and identity["nlevp1"] != want_nlevp1:
        got.append(f"vertres {identity['vertres']:g} (nlevp1 {identity['nlevp1']}), the run needs "
                   f"{float(vertres):g} (nlevp1 {want_nlevp1})")
    if want_nmlat and identity.get("nmlat") and identity["nmlat"] != want_nmlat:
        got.append(f"nres_grid {identity['nres_grid']} (nmlat {identity['nmlat']}), the run needs "
                   f"{int(float(nres_grid))} (nmlat {want_nmlat})")
    if got:
        problems.append(f"{role} {path} was built for " + "; ".join(got)
                        + " - the model stops when its start-up history is on another grid; "
                          "rebuild it for this grid (--compile / --onlycompile)")
    return None, problems


def git_state(repo):
    """{head, short, branch, dirty} of the git checkout at repo; None values when git cannot tell."""
    import subprocess
    out = {"head": None, "short": None, "branch": None, "dirty": None}
    if not repo:
        return out

    def git(*args):
        try:
            p = subprocess.run(["git", "-C", str(repo), "-c", "safe.directory=*", *args],
                               capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError):
            return None
        return p.stdout.strip() if p.returncode == 0 else None

    head = git("rev-parse", "HEAD")
    if not head:
        return out
    out["head"] = head
    out["short"] = git("rev-parse", "--short", "HEAD")
    out["branch"] = git("rev-parse", "--abbrev-ref", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=no")
    out["dirty"] = None if status is None else bool(status)
    return out


def svn_revision_stamp(state):
    """The Make.env SVN_REVISION: git short hash, '+' if dirty; at most 16 characters (the
    history attribute length)."""
    short = (state or {}).get("short")
    if not short:
        return "unknown"
    plus = "+" if state.get("dirty") else ""
    return short[:16 - len(plus)] + plus


# The options dict holds paths for tiegcmrun's own I/O. Generated job files name paths relative to
# the job's workdir and saved JSONs relative to their own directory, inside the run root only.

RUN_DIR_KEYS = ("parentdir", "execdir", "workdir", "histdir")
RUN_FILE_KEYS = ("modelexe", "coupled_modelexe", "log_file", "input_file")


def _path_value(value):
    """True for a string that names a path (not unset, 'gen' or a $-value)."""
    return (isinstance(value, str) and value.strip().lower() not in ("", "none", "null", "gen")
            and "$" not in value)


def run_root(data):
    """parentdir when it contains workdir, else workdir; None without a workdir."""
    if not _path_value(data.get("workdir")):
        return None
    workdir = os.path.abspath(data["workdir"])
    parent = data.get("parentdir")
    if _path_value(parent):
        parent = os.path.abspath(parent)
        if os.path.commonpath([parent, workdir]) == parent:
            return parent
    return workdir


def job_path(path, data, dot=True):
    """path as the job names it: relative to workdir inside the run root, else absolute."""
    root = run_root(data)
    if root is None:
        return path
    return spell(path, data["workdir"], root, dot=dot)


def history_dir(data):
    """histdir as the namelist OUTPUT / SECOUT name it."""
    return job_path(data["histdir"], data)


def job_inp(inp, data):
    """A copy of inp with its file keys as the job names them."""
    out = dict(inp)
    for key in INP_FILE_KEYS:
        if _path_value(out.get(key)):
            out[key] = job_path(out[key], data)
    return out


def json_paths(options, json_dir):
    """A copy of options for saving in json_dir: run paths inside the run root become relative
    to json_dir."""
    out = copy.deepcopy(options)
    data = (out.get("model") or {}).get("data")
    if not isinstance(data, dict):
        return out
    root = run_root(data)
    if root is None:
        return out
    rel = lambda p: spell(p, json_dir, root)        # noqa: E731
    inp = out.get("inp")
    if isinstance(inp, dict):
        for key in INP_FILE_KEYS:
            if _path_value(inp.get(key)):
                inp[key] = rel(inp[key])
    for key in RUN_DIR_KEYS + RUN_FILE_KEYS:
        if _path_value(data.get(key)):
            data[key] = rel(data[key])
    return out


def resolve_json_paths(options, json_dir):
    """Resolve the relative run paths of a loaded options JSON against json_dir, in place."""
    data = (options.get("model") or {}).get("data")
    base = os.path.abspath(json_dir)
    fix = lambda p: (os.path.normpath(os.path.join(base, p))          # noqa: E731
                     if _path_value(p) and not os.path.isabs(p) else p)
    if isinstance(data, dict):
        for key in RUN_DIR_KEYS + RUN_FILE_KEYS:
            if key in data:
                data[key] = fix(data[key])
    inp = options.get("inp")
    if isinstance(inp, dict):
        for key in INP_FILE_KEYS:
            if key in inp:
                inp[key] = fix(inp[key])
    return options


OPTION_DESCRIPTIONS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "options_description.json")
NONE_SPELLINGS = ("", "none", "null", "[none]", "[null]")
_NAME_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")
_HMS_RE = re.compile(r"^(\d+):(\d{2}):(\d{2})$")
_TRUE = ("true", "t", "yes", "y", "1", ".true.", ".t.")
_FALSE = ("false", "f", "no", "n", "0", ".false.", ".f.")


class ValidsError(ValueError):
    """A value of the right type that is not one of the field's valids."""


def is_unset(value):
    """True for None, a None spelling, or an empty / [None] list."""
    if value is None:
        return True
    if isinstance(value, (list, tuple)):
        return all(is_unset(v) for v in value)
    return isinstance(value, str) and value.strip().lower() in NONE_SPELLINGS


@lru_cache(maxsize=None)
def _field_index():
    with open(OPTION_DESCRIPTIONS_FILE, encoding="utf-8") as f:
        od = json.load(f)
    out = []

    def walk(node, prefix):
        for name, meta in node.items():
            if not isinstance(meta, dict):
                continue
            path = f"{prefix}.{name}" if prefix else name
            if "LEVEL" in meta:
                out.append((path, name, meta))
            else:
                walk(meta, path)
    walk(od, "")
    return tuple(out)


def field_meta(prefix=None):
    """{leaf name: description} of options_description.json under the dotted prefix."""
    out = {}
    for path, name, meta in _field_index():
        if prefix is None or path.startswith(prefix + "."):
            out.setdefault(name, meta)
    return out


def keys_of_type(kind, prefix=None):
    kinds = (kind,) if isinstance(kind, str) else tuple(kind)
    return tuple(n for n, m in field_meta(prefix).items() if m.get("type") in kinds)


INP_FILE_KEYS = keys_of_type("file", "inp")


def anchor_path(value, run_dir):
    """A relative path resolved against run_dir; other values unchanged."""
    if not run_dir or not _path_value(value):
        return value
    path = os.path.expanduser(value.strip())
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(run_dir, path))


def _number(name, value):
    """float(value) for a finite number, Fortran 1.5D-5 included."""
    if isinstance(value, bool):
        raise ValueError(f"{name} = {value!r} is not a number")
    try:
        x = float(value) if not isinstance(value, str) else float(value.strip().lower().replace("d", "e"))
    except (TypeError, ValueError):
        raise ValueError(f"{name} = {value!r} is not a number") from None
    if x != x or x in (float("inf"), float("-inf")):
        raise ValueError(f"{name} = {value!r} is not a finite number")
    return x


def _integer(name, value):
    if isinstance(value, bool):
        raise ValueError(f"{name} = {value!r} is not a whole number")
    if isinstance(value, int):
        return value
    try:
        x = _number(name, value)
    except ValueError:
        raise ValueError(f"{name} = {value!r} is not a whole number") from None
    if not x.is_integer():
        raise ValueError(f"{name} = {value!r} is not a whole number")
    return int(x)


def _boolean(name, value):
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    if s in _TRUE:
        return True
    if s in _FALSE:
        return False
    raise ValueError(f"{name} = {value!r} is not true or false")


def _bounds(name, x, desc, shown):
    lo, hi = desc.get("min"), desc.get("max")
    if lo is not None:
        if desc.get("exclusive_min") and not x > lo:
            raise ValueError(f"{name} = {shown!r} must be > {lo}")
        if not desc.get("exclusive_min") and x < lo:
            raise ValueError(f"{name} = {shown!r} must be >= {lo}")
    if hi is not None and x > hi:
        raise ValueError(f"{name} = {shown!r} must be <= {hi}")


def _check_valids(name, value, desc):
    valids = desc.get("valids")
    if not valids:
        return
    if str(value).strip() in [str(v).strip() for v in valids]:
        return
    if desc.get("type") in ("int", "float"):
        try:
            x = float(str(value).strip().lower().replace("d", "e"))
            if any(abs(x - float(v)) < 1e-9 for v in valids):
                return
        except (TypeError, ValueError):
            pass
    raise ValidsError(f"{name} = {value!r} is not one of {' | '.join(map(str, valids))}")


def _dhms(name, value):
    parts = value if isinstance(value, (list, tuple)) else str(value).replace("'", "").replace(",", " ").split()
    try:
        ints = [_integer(name, p) for p in parts]
    except ValueError:
        ints = None
    if ints is None or len(ints) != 4:
        raise ValueError(f"{name} = {value!r} is not 4 integers D H M S")
    if not any(ints):
        raise ValueError(f"{name} = {value!r} must not be all zero")
    return " ".join(map(str, ints))


def _file_like(name, value, kind, run_dir, check_exists=True):
    s = str(value).strip()
    if s.lower() == "gen" or "$" in s:
        return s
    path = anchor_path(s, run_dir) if run_dir else os.path.expanduser(s)
    if kind == "path" or not check_exists:
        return path
    if kind == "dir":
        if os.path.isdir(path):
            return path
        raise ValueError(f"{name} = {value!r} is not a directory")
    if os.path.isfile(path):
        return path
    if os.path.isdir(path):
        raise ValueError(f"{name} = {value!r} is a directory, not a file")
    if os.sep not in s:
        found = find_file(s, tiegcm_env("TIEGCMDATA") or "")
        if found:
            return found
    raise ValueError(f"{name} = {value!r} not found")


def coerce(name, value, desc, run_dir=None, check_exists=True):
    """The canonical value for field desc of options_description.json; unset values give None.

    Raises ValueError (ValidsError for a valids miss) naming the key. Floats are checked but kept
    as given; relative file paths are anchored to run_dir."""
    kind = (desc or {}).get("type")
    if kind is None:
        return value
    if is_unset(value):
        return None
    key = name.rsplit(".", 1)[-1]
    if kind == "int":
        out = _integer(name, value)
        _bounds(name, out, desc, value)
    elif kind == "float":
        _bounds(name, _number(name, value), desc, value)
        out = value.strip() if isinstance(value, str) else value
    elif kind == "bool":
        return _boolean(name, value)
    elif kind in ("str", "name"):
        out = str(value).strip()
        if kind == "name" and not _NAME_RE.match(out):
            raise ValueError(f"{name} = {value!r}: use only letters, digits and _ . + -")
    elif kind == "iso":
        from namelist_solver import parse_run_datetime, RUN_DATETIME_FORMAT
        try:
            out = parse_run_datetime(value).strftime(RUN_DATETIME_FORMAT)
        except ValueError:
            raise ValueError(f"{name} = {value!r} is not yyyy-mm-ddThh:mm:ss") from None
    elif kind == "dhms":
        return _dhms(name, value)
    elif kind == "hms":
        m = _HMS_RE.match(str(value).strip())
        if not m or int(m.group(2)) > 59 or int(m.group(3)) > 59:
            raise ValueError(f"{name} = {value!r} is not HH:MM:SS")
        return str(value).strip()
    elif kind in ("file", "dir", "path"):
        return _file_like(name, value, kind, run_dir, check_exists)
    elif kind == "lines":
        return [x for x in as_list(value)] or [None]
    elif kind.endswith("_list"):
        if key in NAMELIST_ARRAYS:
            try:
                out = namelist_array(key, value)
            except ValueError as e:
                raise ValueError(f"{name.rsplit('.', 1)[0] + '.' if '.' in name else ''}{e}") from None
        else:
            out = [str(x).strip("'\"") for x in as_list(value, tokens=True) if x is not None]
        if kind in ("int_list", "float_list"):
            for x in out or []:
                _bounds(name, _number(name, x), desc, x)
        return out
    else:
        return value
    _check_valids(name, out, desc)
    return out
