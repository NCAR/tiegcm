"""Arrow-navigable option wizard for tiegcmrun.

Asks the same questions as the linear prompts (tiegcmrun.get_run_option) and returns the same
options dict.
"""
import copy
import datetime
import json
import os
import re
import socket
import warnings
from fractions import Fraction

LEVELS = {"BENCH": -1, "BASIC": 0, "INTERMEDIATE": 1, "EXPERT": 2}


def tiegcm_env(name):
    from misc import tiegcm_env as env
    return env(name)


def _misc():
    import misc
    return misc


class _FormUnavailable(Exception):
    """The wizard cannot run (no prompt_toolkit or no TTY); the caller falls back."""


class _FormCancelled(_FormUnavailable):
    """Ctrl-C: the caller exits without writing and does not fall back. Catch it first."""


class _KeepDefault(Exception):
    """Raised by a derive to keep the field's static default (not a failure)."""


class Field:
    __slots__ = ("name", "path", "level", "default", "valids", "prompt", "section",
                 "description", "warning", "meta")

    def __init__(self, name, path, level, default, valids, prompt, section,
                 description=None, warning=None, meta=None):
        self.name, self.path, self.level = name, path, level
        self.default, self.valids, self.prompt, self.section = default, valids, prompt, section
        self.description, self.warning = description, warning
        self.meta = meta or {}                  # the raw options_description.json entry


class FormState:
    """Field values keyed by path; re-derives the defaults the user has not overridden.

    derive: {path: (fn(values) -> value, deps)}; valid: {path: fn(values, value) -> error or None}.
    """
    def __init__(self, fields, derive, valid, mode="EXPERT", context=None, cond=None, warn=None):
        self.fields = fields
        self.by_path = {f.path: f for f in fields}
        self.derive = derive
        self.valid = valid
        self.mode = mode
        self.cond = cond                        # cond(field, values) -> shown?
        self.warn = warn or {}                  # path -> fn(values) -> warning or None
        self.context = context or {}            # __engage/__benchmark/__coupling flags
        self.values = {f.path: ("" if f.default is None else f.default) for f in fields}
        self.overridden = set()
        self.derive_errors = {}
        self.rederive()

    def _all(self):
        d = dict(self.values)
        d.update(self.context)
        return d

    def visible(self):
        allv = self._all()
        # a blank parentdir asks for the run dirs even in BASIC
        pd_blank = (not _is_set(allv.get("model.data.parentdir"))) and (not allv.get("__engage"))
        out = []
        for f in self.fields:
            if pd_blank and f.path in ("model.data.execdir", "model.data.workdir",
                                       "model.data.histdir"):
                out.append(f)
                continue
            if LEVELS.get(f.level, 2) > LEVELS[self.mode]:
                continue
            if self.cond is not None and not self.cond(f, allv):
                continue
            out.append(f)
        return out

    def set(self, path, value, normalize=True):
        if normalize:
            value = self._normalize(path, value)
            # a cleared parentdir means the default, not the flat layout
            if path == "model.data.parentdir" and value == "" and path in self.derive:
                self.reset(path)
                return
        self.values[path] = value
        self.overridden.add(path)
        self.rederive()

    def reset(self, path):
        """Ctrl-R: restore the derived or static default."""
        self.overridden.discard(path)
        f = self.by_path.get(path)
        if f is not None:
            self.values[path] = _static_default(f)
        self.rederive()

    def _normalize(self, path, value):
        """Clean up a typed answer as tiegcmrun.get_run_option does."""
        f = self.by_path.get(path)
        if f is None or not isinstance(value, str) or f.name in _LIST_FIELDS:
            return value                          # parsed by misc.as_list on emit
        s = value.strip()
        if f.name in _FOURVAR_FIELDS:
            if s.lower() == "none":
                # only the segment length is optional; elsewhere the validator rejects 'none'
                return [None] if f.name == "segment" else s
            try:                                  # '1,0,0,0' / "'0 6 0 0'" -> '1 0 0 0'
                ints = [int(x) for x in s.replace("'", "").replace(",", " ").split()]
            except ValueError:
                return s
            return " ".join(map(str, ints)) if ints else s
        if s.lower() == "none":
            return None
        if f.name == "vertres":                   # accepted as 1/2|1/4|1/8|1/16
            try:
                return float(Fraction(s))
            except (ValueError, ZeroDivisionError):
                return s
        if (f.meta.get("type") in ("file", "dir", "path") and f.name not in _misc().RUN_DIR_KEYS
                and f.section != "simulation"):
            # relative to the run folder (parentdir, else workdir)
            root = _misc().run_root({"parentdir": self.values.get("model.data.parentdir"),
                                     "workdir": self.values.get("model.data.workdir")})
            s = _misc().anchor_path(s, os.path.abspath(root)) if root else s
        if (f.name in _FILE_FIELDS and s and s.lower() != "gen"
                and not os.path.isfile(s) and not os.path.isdir(s)
                and os.sep not in value.strip()):
            from misc import find_file
            found = find_file(value.strip(), tiegcm_env("TIEGCMDATA") or "")
            if found:
                return str(found)
        return s

    def rederive(self):
        # prompt order puts dependencies first, so one forward pass suffices
        self.derive_errors = {}
        for f in self.fields:
            if f.path in self.overridden:
                continue
            spec = self.derive.get(f.path)
            if not spec:
                dflt = _static_default(f)
                if self.values[f.path] != dflt:
                    self.values[f.path] = dflt
                continue
            fn, _deps = spec
            try:
                self.values[f.path] = fn(self._all())
            except _KeepDefault:
                self.values[f.path] = _static_default(f)
            except Exception as e:
                self.derive_errors[f.path] = str(e)   # keep the prior value

    def errors(self):
        out = {}
        for f in self.visible():
            vfn = self.valid.get(f.path)
            if vfn:
                e = vfn(self._all(), self.values[f.path])
                if e:
                    out[f.path] = e
                    continue
            # a failed derive leaves a stale value
            if f.path in self.derive_errors and f.path not in self.overridden:
                out[f.path] = f"could not derive: {self.derive_errors[f.path]}"
                continue
            v = str(self.values[f.path]).strip()
            listed = str(f.meta.get("type", "")).endswith("_list")   # valids of each element
            if f.valids and v != "" and not listed:
                if _is_bool_field(f):
                    if _coerce_bool(v) is None:
                        out[f.path] = "must be true or false"
                    continue
                vv = [str(x).strip() for x in f.valids]
                ok = v in vv
                if not ok:                        # "2.50" matches 2.5
                    fv = _try_float(v)
                    ok = fv is not None and fv in [g for g in (_try_float(x) for x in vv)
                                                   if g is not None]
                if not ok:
                    out[f.path] = f"must be one of {f.valids}"
                    continue
            if v != "":
                from misc import coerce
                try:
                    coerce(f.name, self.values[f.path], f.meta)
                except ValueError as e:
                    out[f.path] = str(e)
        return out

    def warnings_for(self, path, live=None):
        """Return the value-dependent warning for `path`, or None.

        `live` is the uncommitted input text, used instead of the committed value.
        """
        fn = self.warn.get(path)
        if not fn:
            return None
        allv = self._all()
        if live is not None:
            allv[path] = live
        try:
            return fn(allv) or None
        except Exception:
            return None

    def options_dict(self):
        """Return the values as the nested options dict.

        Bool fields become Python bools: the template's `is true` tests would drop a 'true' string.
        """
        out = {}
        allv = self._all()
        for f in self.fields:
            cur = out
            parts = f.path.split(".")
            for p in parts[:-1]:
                cur = cur.setdefault(p, {})
            val = self.values[f.path]
            # a hidden forcing field is omitted even if it was set before it was hidden
            if (f.name in _FORCING_SKIP + _GSWM_NM_FIELDS and self.cond is not None
                    and not self.cond(f, allv)):
                val = None
            elif _emits_bool(f):
                cb = _coerce_bool(val)
                if cb is not None:
                    val = cb
            # List fields must be lists (the templates iterate them); a cleared box is [None].
            if isinstance(val, str) and f.name in _LIST_FIELDS and not (
                    val.strip() == "" and f.default is None):
                from misc import as_list, ascii_quotes, secflds_list
                if f.name == "SECFLDS":
                    val = as_list(ascii_quotes(val), tokens=True)
                    try:
                        val = secflds_list(val)[0]
                    except ValueError:
                        pass
                else:
                    val = as_list(val, tokens=f.name in _TOKEN_LIST_FIELDS)
            # unset is None: the template would write an empty namelist line for ""
            elif isinstance(val, str) and val.strip() == "":
                val = None
            cur[parts[-1]] = val
        return out


def _static_default(f):
    """Return a copy of the field's static default ('' for null)."""
    return "" if f.default is None else copy.deepcopy(f.default)


def _try_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


_BOOL_TRUE = {"yes", "y", "true", "t", "1"}
_BOOL_FALSE = {"no", "n", "false", "f", "0"}

# Token lists split on commas/spaces too; the other list fields hold one command per line or ';'.
_ARRAY_FIELDS = tuple(_misc().NAMELIST_ARRAYS)
_LIST_FIELDS = ("SECFLDS", "other_input", "other_job", "local_modules", "job_chain") + _ARRAY_FIELDS
_TOKEN_LIST_FIELDS = ("SECFLDS",) + _ARRAY_FIELDS

# warned about against the uncommitted input text
_LIVE_WARN = ("model.data.modelexe", "model.data.coupled_modelexe")

_FOURVAR_FIELDS = _misc().keys_of_type("dhms")
_FILE_FIELDS = _misc().INP_FILE_KEYS

# the model refuses non-migrating GSWM files at 5 degrees
_GSWM_NM_FIELDS = ("GSWM_NM_DI_NCFILE", "GSWM_NM_SDI_NCFILE")

# omitted when GPI/IMF/KP/POTENTIAL_MODEL make them unused
_FORCING_SKIP = ("KP", "POWER", "CTPOTEN", "F107", "F107A", "IMF_NCFILE",
                 "BXIMF", "BYIMF", "BZIMF", "SWDEN", "SWVEL")


def _is_bool_field(f):
    """True for a scalar [true,false] field (logical arrays are list fields)."""
    return (bool(f.valids) and {str(x).strip().lower() for x in f.valids} == {"true", "false"}
            and f.meta.get("type", "bool") == "bool")


def _coerce_bool(value):
    s = str(value).strip().lower()
    if s in _BOOL_TRUE:
        return True
    if s in _BOOL_FALSE:
        return False
    return None


def _emits_bool(f):
    """True for a [true,false] field; [0,1] flags are Fortran integers and must stay 1/0."""
    return _is_bool_field(f)


def _iso(s):
    from namelist_solver import parse_run_datetime
    return parse_run_datetime(s)


def load_tiegcmrun_fields(option_descriptions):
    """Flatten options_description.json into a Field list in prompt (dependency) order."""
    fields = []

    def walk(node, prefix, section):
        for name, meta in node.items():
            if isinstance(meta, dict) and ("LEVEL" in meta and ("prompt" in meta or "default" in meta)):
                path = f"{prefix}.{name}" if prefix else name
                fields.append(Field(name, path, meta.get("LEVEL", "EXPERT"),
                                     meta.get("default"), meta.get("valids"),
                                     meta.get("prompt", name), section,
                                     meta.get("description"), meta.get("warning"), meta))
            elif isinstance(meta, dict):
                walk(meta, f"{prefix}.{name}" if prefix else name, section or name)

    for top in ("simulation", "model", "inp", "job"):
        if top in option_descriptions:
            walk(option_descriptions[top], top, top)
    return fields


def _derive_registry():
    """Return {path: (fn(values) -> default, deps)}, built on tiegcmrun's own functions."""
    from misc import (resolution_solver, get_mtime, select_source_defaults, find_file,
                      select_latest_gpi, file_covers, mres_to_nres_grid, default_make_file,
                      gswm_allowed, he_coefs_pattern, GSWM_PATTERNS, history_dir)
    from namelist_solver import (inp_pri_date, inp_prihist, inp_sechist, inp_mxhist,
                                  inp_pri_out, inp_sec_out, inp_sec_date, cadence_step,
                                  lbc_other_set)
    from output_solver import queue_walltime_default
    from engage_solver import COUPLED_CAPS

    def res(i):
        return lambda v: resolution_solver(v["model.specification.horires"])[i]

    def pri_date(i):
        return lambda v: " ".join(map(str, inp_pri_date(v["inp.start_time"], v["inp.stop_time"])[i])) \
            if i >= 2 else inp_pri_date(v["inp.start_time"], v["inp.stop_time"])[i]

    def _cadence_step(v, *steps):
        return cadence_step(v.get("__engage_cfg"), *steps)

    def _seg(v):
        s = str(v["inp.segment"]).strip()
        if s in ("", "[none]", "[None]", "[null]", "none", "None", "[]"):
            return None
        try:
            seg = [int(x) for x in s.split()]
        except ValueError:
            return None
        return seg if any(seg) else None

    def _source(v):
        if v.get("__benchmark"):
            return find_file(f'{v["__benchmark_name"]}.nc', tiegcm_env("TIEGCMDATA") or "")
        return select_source_defaults(
            {"inp": {"start_time": v["inp.start_time"],
                     "solar_flux_level": v["inp.solar_flux_level"]}}, {})

    _mtime_cache = {}

    def _src_start(v):
        # first SOURCE history at the run start's time of day (input.F requires it);
        # memoized because get_mtime opens a netCDF
        src = str(v.get("inp.SOURCE", "")).strip()
        if src in ("", "None", "none"):
            return ""
        if src not in _mtime_cache:
            with warnings.catch_warnings():       # stderr would corrupt the TUI
                warnings.simplefilter("ignore")
                _mtime_cache[src] = get_mtime(src)
        mtimes = _mtime_cache[src]
        if not v.get("__benchmark"):
            try:
                st = _iso(v["inp.start_time"])
                mtimes = [m for m in mtimes if list(m[1:4]) == [st.hour, st.minute, st.second]] or mtimes
            except (KeyError, TypeError, ValueError):
                pass
        return " ".join(map(str, mtimes[0]))

    def _run_name(v):
        return (f'{v["simulation.job_name"]}_{v["model.specification.horires"]}'
                f'x{v["model.specification.vertres"]}')

    def _hpc_default(v):
        eng = v.get("__engage_cfg")
        if eng and eng.get("hpc_system"):
            return eng["hpc_system"]
        from output_solver import CONFIG_DIR
        from jobgen import JobGen
        return JobGen(CONFIG_DIR, hostname=socket.gethostname()).name

    def _eng_or_keep(key):
        """Take `key` from the engage payload when coupled, else keep the static default."""
        def fn(v):
            cfg = v.get("__engage_cfg")
            if isinstance(cfg, dict) and cfg.get(key) is not None:
                return cfg[key]
            raise _KeepDefault(key)
        return fn

    def _account(machine):
        # engage's account, else the machine default, else the user's only project code
        keep = _eng_or_keep("project_code")

        def fn(v):
            try:
                return keep(v)
            except _KeepDefault:
                from output_solver import CONFIG_DIR
                from jobgen import JobGen, user_projects
                m = JobGen(CONFIG_DIR, machine=machine).m
                if m.get("account_groups") and m.get("account_default") is None:
                    choices = user_projects()
                    if len(choices) == 1:
                        return choices[0]
                raise
        return fn

    def _walltime(machine):
        def fn(v):
            cfg = v.get("__engage_cfg")
            if isinstance(cfg, dict) and cfg.get("walltime") is not None:
                return cfg["walltime"]
            queue = str(v.get(f"job.{machine}.queue") or "").strip()
            wall = queue_walltime_default(machine, queue) if queue else None
            if wall is None:
                raise _KeepDefault("walltime")
            return wall
        return fn

    def _coupled_cap(key):
        def fn(v):
            if v.get("__engage"):
                return COUPLED_CAPS[key]
            raise _KeepDefault(key)
        return fn

    def _make(v):
        return default_make_file(v["model.data.modeldir"], v["simulation.hpc_system"])

    _gpi_cache = {}

    def _gpi(v):
        # newest local GPI if it covers the window, else "gen" (generated on accept)
        start, stop = v["inp.start_time"], v["inp.stop_time"]
        if "latest" not in _gpi_cache:
            _gpi_cache["latest"] = select_latest_gpi(tiegcm_env("TIEGCMDATA") or "")
        latest = _gpi_cache["latest"]
        if not latest:
            return "gen"
        key = (latest, start, stop)
        if key not in _gpi_cache:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _gpi_cache[key] = file_covers(latest, start, stop)
        return latest if _gpi_cache[key] else "gen"

    def _gswm(key, pattern):
        # float(): files are named by the float resolution ('..._5.0d_99km.nc')
        def fn(v):
            horires = float(v["model.specification.horires"])
            if not gswm_allowed(key, horires) or lbc_other_set(_inp(v)):
                return ""
            return find_file(pattern.format(horires=horires), tiegcm_env("TIEGCMDATA") or "") or ""
        return fn

    def _he_coefs(v):
        found = find_file(he_coefs_pattern(v["model.specification.horires"]),
                          tiegcm_env("TIEGCMDATA") or "")
        return found if found is not None else ""

    def _rundir(sub):
        # <parentdir>/<sub>; no parentdir: one flat dir; coupled: everything in parentdir
        def fn(v):
            pd = v["model.data.parentdir"]
            if v.get("__engage"):
                return pd
            if not _is_set(pd):
                return "." if sub == "exec" else v["model.data.execdir"]
            return os.path.join(str(pd), sub)
        return fn

    def _histdir_job(v):
        return history_dir({"parentdir": v.get("model.data.parentdir"),
                            "workdir": v["model.data.workdir"], "histdir": v["model.data.histdir"]})

    def _res(machine, index):
        def fn(v):
            from misc import select_resource_defaults
            probe = {"simulation": {"hpc_system": machine},
                     "model": {"specification": {"horires": v["model.specification.horires"],
                                                 "nres_grid": v.get("model.specification.nres_grid")}},
                     "job": {"resource": {"model": str(v.get(f"job.{machine}.resource.model")).strip()},
                             "queue": str(v.get(f"job.{machine}.queue") or "").strip()}}
            try:
                out = select_resource_defaults(probe, {})[index]
            except (ValueError, TypeError, KeyError):
                return (1, 128, 128)[index]         # an aitken model the table does not know
            if index == 2:                          # no more ranks per node than cores
                try:
                    out = min(int(out), int(v.get(f"job.{machine}.resource.ncpus")))
                except (TypeError, ValueError):
                    pass
            return out
        return fn

    R = {
        "simulation.job_name": (lambda v: (v.get("__engage_cfg") or {}).get("job_name")
                                or v.get("__benchmark_name") or "tiegcm", []),
        "model.specification.vertres":   (res(0), ["model.specification.horires"]),
        "model.specification.mres":      (res(1), ["model.specification.horires"]),
        "model.specification.nres_grid": (lambda v: mres_to_nres_grid(v["model.specification.mres"]),
                                          ["model.specification.mres"]),
        "inp.STEP":                      (res(3), ["model.specification.horires"]),
        "model.data.modeldir":  (lambda v: tiegcm_env("TIEGCMHOME") or "", []),
        "model.data.parentdir": (lambda v: (v.get("__engage_cfg") or {}).get("parentdir")
                                  or ".", []),
        "model.data.execdir":   (_rundir("exec"), ["model.data.parentdir"]),
        "model.data.workdir":   (_rundir("stdout"), ["model.data.parentdir", "model.data.execdir"]),
        "model.data.histdir":   (_rundir("hist"), ["model.data.parentdir", "model.data.execdir"]),
        "model.data.tgcmdata":  (lambda v: tiegcm_env("TIEGCMDATA") or "", []),
        "model.data.log_file":  (lambda v: f'{v["model.data.workdir"]}/{v["simulation.job_name"]}.out',
                                 ["model.data.workdir", "simulation.job_name"]),
        "model.data.make":      (_make, ["model.data.modeldir", "simulation.hpc_system"]),
        "model.data.modelexe":  (lambda v: os.path.join(v["model.data.execdir"], "tiegcm.exe"),
                                 ["model.data.execdir"]),
        "model.data.coupled_modelexe": (lambda v: os.path.join(v["model.data.execdir"], "tiegcm.x"),
                                 ["model.data.execdir"]),
        "inp.LABEL": (lambda v: f'{v["simulation.job_name"]}_{v["model.specification.horires"]}'
                      f'x{v["model.specification.vertres"]}',
                      ["simulation.job_name", "model.specification.horires",
                       "model.specification.vertres"]),
        "inp.SOURCE":       (_source, ["inp.start_time", "inp.solar_flux_level"]),
        "inp.SOURCE_START": (_src_start, ["inp.SOURCE", "inp.start_time"]),
        "inp.START_YEAR": (pri_date(0), ["inp.start_time", "inp.stop_time"]),
        "inp.START_DAY":  (pri_date(1), ["inp.start_time", "inp.stop_time"]),
        "inp.PRISTART":   (pri_date(2), ["inp.start_time", "inp.stop_time"]),
        "inp.PRISTOP":    (pri_date(3), ["inp.start_time", "inp.stop_time"]),
        "inp.secondary_start_time": (lambda v: v["inp.start_time"], ["inp.start_time"]),
        "inp.secondary_stop_time":  (lambda v: v["inp.stop_time"], ["inp.stop_time"]),
        # cadences are multiples of STEP; PRIHIST is asked before STEP, so it uses the default
        "inp.PRIHIST": (lambda v: " ".join(map(str, inp_prihist(
            [int(x) for x in v["inp.PRISTART"].split()], [int(x) for x in v["inp.PRISTOP"].split()], _seg(v),
            _cadence_step(v, resolution_solver(v["model.specification.horires"])[3])))),
            ["inp.PRISTART", "inp.PRISTOP", "inp.segment", "model.specification.horires"]),
        "inp.SECHIST": (lambda v: " ".join(map(str, inp_sechist(
            *inp_pri_date(v["inp.secondary_start_time"], v["inp.secondary_stop_time"])[2:4], _seg(v),
            _cadence_step(v, v.get("inp.STEP"), resolution_solver(v["model.specification.horires"])[3])))),
            ["inp.secondary_start_time", "inp.secondary_stop_time", "inp.segment", "inp.STEP",
             "model.specification.horires"]),
        "inp.MXHIST_PRIM": (lambda v: inp_mxhist(v["inp.start_time"], v["inp.stop_time"],
            [int(x) for x in v["inp.PRIHIST"].split()], None, _seg(v))[0],
            ["inp.PRIHIST", "inp.start_time", "inp.stop_time", "inp.segment"]),
        "inp.MXHIST_SECH": (lambda v: inp_mxhist(v["inp.start_time"], v["inp.stop_time"],
            [int(x) for x in v["inp.SECHIST"].split()], None, _seg(v))[0],
            ["inp.SECHIST", "inp.start_time", "inp.stop_time", "inp.segment"]),
        "inp.SECSTART": (lambda v: " ".join(map(str, inp_sec_date(
            v["inp.secondary_start_time"], v["inp.secondary_stop_time"], [int(x) for x in v["inp.SECHIST"].split()])[0])),
            ["inp.secondary_start_time", "inp.secondary_stop_time", "inp.SECHIST"]),
        "inp.SECSTOP": (lambda v: " ".join(map(str, inp_sec_date(
            v["inp.secondary_start_time"], v["inp.secondary_stop_time"], [int(x) for x in v["inp.SECHIST"].split()])[1])),
            ["inp.secondary_start_time", "inp.secondary_stop_time", "inp.SECHIST"]),
        "simulation.hpc_system": (_hpc_default, []),
        "model.specification.segmentation": (
            lambda v: bool(_is_set(v["inp.segment"]) and not _is_set(v["model.data.input_file"])),
            ["inp.segment", "model.data.input_file"]),
        "inp.OUTPUT": (lambda v: inp_pri_out(
            v["inp.start_time"], v["inp.stop_time"], [int(x) for x in str(v["inp.PRIHIST"]).split()],
            int(v["inp.MXHIST_PRIM"]), 0, _histdir_job(v), _run_name(v))[0],
            ["inp.start_time", "inp.stop_time", "inp.PRIHIST", "inp.MXHIST_PRIM",
             "model.data.histdir", "model.data.workdir", "model.data.parentdir",
             "simulation.job_name", "model.specification.horires", "model.specification.vertres"]),
        "inp.SECOUT": (lambda v: inp_sec_out(
            v["inp.secondary_start_time"], v["inp.secondary_stop_time"],
            [int(x) for x in str(v["inp.SECHIST"]).split()], int(v["inp.MXHIST_SECH"]), 0,
            _histdir_job(v), _run_name(v))[0],
            ["inp.secondary_start_time", "inp.secondary_stop_time", "inp.SECHIST",
             "inp.MXHIST_SECH", "model.data.histdir", "model.data.workdir", "model.data.parentdir",
             "simulation.job_name", "model.specification.horires", "model.specification.vertres"]),
        "inp.GPI_NCFILE": (_gpi, ["inp.start_time", "inp.stop_time"]),
        **{"inp." + key: (_gswm(key, pattern), ["model.specification.horires"])
           for key, pattern in GSWM_PATTERNS.items()},
        "inp.HE_COEFS_NCFILE": (_he_coefs, ["model.specification.horires"]),
        **{"inp." + key: (_coupled_cap(key), []) for key in COUPLED_CAPS},
        "job.aitken.resource.model": (_eng_or_keep("model"), []),
        "job.derecho.project_code": (_account("derecho"), []),
        "job.derecho.queue":        (_eng_or_keep("queue"), []),
        "job.derecho.job_priority": (_eng_or_keep("job_priority"), []),
        "job.derecho.walltime":     (_walltime("derecho"), ["job.derecho.queue"]),
        "job.aitken.project_code":  (_account("aitken"), []),
        "job.aitken.queue":         (_eng_or_keep("queue"), []),
        "job.aitken.walltime":      (_walltime("aitken"), ["job.aitken.queue"]),
        "job.aitken.group_list":    (_eng_or_keep("group_list"), []),
        "job.derecho.resource.select":   (_res("derecho", 0), ["model.specification.horires", "job.derecho.queue"]),
        "job.derecho.resource.ncpus":    (_res("derecho", 1), ["model.specification.horires", "job.derecho.queue"]),
        "job.derecho.resource.mpiprocs": (_res("derecho", 2), ["model.specification.horires", "job.derecho.queue",
                                                               "job.derecho.resource.ncpus"]),
        "job.derecho.nprocs": (lambda v: int(v["job.derecho.resource.select"]) * int(v["job.derecho.resource.mpiprocs"]),
            ["job.derecho.resource.select", "job.derecho.resource.mpiprocs"]),
        "job.aitken.resource.select":   (_res("aitken", 0), ["model.specification.horires", "job.aitken.resource.model", "job.aitken.queue"]),
        "job.aitken.resource.ncpus":    (_res("aitken", 1), ["model.specification.horires", "job.aitken.resource.model", "job.aitken.queue"]),
        "job.aitken.resource.mpiprocs": (_res("aitken", 2), ["model.specification.horires", "job.aitken.resource.model", "job.aitken.queue",
                                                             "job.aitken.resource.ncpus"]),
        "job.aitken.nprocs": (lambda v: int(v["job.aitken.resource.select"]) * int(v["job.aitken.resource.mpiprocs"]),
            ["job.aitken.resource.select", "job.aitken.resource.mpiprocs"]),
    }
    # any other scheduler machine in config/machines.yaml
    from output_solver import CONFIG_DIR
    from jobgen import JobGen
    machines = JobGen(CONFIG_DIR, machine="derecho").machines
    for name, m in machines.items():
        if name in ("derecho", "aitken") or m.get("scheduler", "none") == "none":
            continue
        pre = f"job.{name}."
        R[pre + "resource.select"] = (_res(name, 0), ["model.specification.horires", pre + "queue"])
        R[pre + "resource.ncpus"] = (_res(name, 1), ["model.specification.horires", pre + "queue"])
        R[pre + "resource.mpiprocs"] = (_res(name, 2), ["model.specification.horires", pre + "queue",
                                                        pre + "resource.ncpus"])
        R[pre + "walltime"] = (_walltime(name), [pre + "queue"])
        R[pre + "nprocs"] = ((lambda v, p=pre: int(v[p + "resource.select"]) * int(v[p + "resource.mpiprocs"])),
                             [pre + "resource.select", pre + "resource.mpiprocs"])
    return R


def _validate_registry():
    """Return {path: fn(values, value) -> error or None}."""
    def hms_match(v, val):
        try:
            ss = [int(x) for x in str(val).split()]
            st = _iso(v["inp.start_time"])
            if (ss[1], ss[2], ss[3]) != (st.hour, st.minute, st.second):
                return (f"time-of-day {ss[1]}:{ss[2]}:{ss[3]} must match run start "
                        f"{st.hour}:{st.minute}:{st.second} (TIE-GCM requires it)")
        except (ValueError, IndexError):
            return None
        return None

    def iso(v, val):
        try:
            _iso(val); return None
        except ValueError:
            return "must be YYYY-MM-DDThh:mm:ss"

    def stop_after(v, val):
        err = iso(v, val)
        if err:
            return err
        try:
            return "must be after start_time" if _iso(val) <= _iso(v["inp.start_time"]) else None
        except (ValueError, KeyError):
            return None

    def seg_required(v, val):
        # input.F caps the day of year at 366/367 for the start year
        from namelist_solver import year_boundary
        try:
            crossing = year_boundary(v["inp.start_time"], v["inp.stop_time"])
            if crossing and not _is_set(val):
                return (f"segmentation required: run crosses the year boundary "
                        f"(stop day {crossing[1]} > {crossing[2]})")
        except Exception:
            return None
        return None

    def sec_window(key):
        from namelist_solver import date_order_problem
        def check(v, val):
            if not _is_set(val):
                return None
            if key == "secondary_start_time":
                return date_order_problem(key, val, v.get("inp.start_time"), "start_time",
                                          strictly_after=False, before=v.get("inp.stop_time"),
                                          before_label="stop_time", strictly_before=True)
            return date_order_problem(key, val, v.get("inp.secondary_start_time"),
                                      "secondary_start_time", before=v.get("inp.stop_time"),
                                      before_label="stop_time")
        return check

    def fourvar(v, val):                            # exactly 4 ints, not all-zero
        s = str(val).strip()
        if isinstance(val, str) and s.lower() == "none":
            return "required: give 4 integers (D H M S)"
        if s == "" or s.lower() in ("none", "[none]", "[null]"):
            return None
        try:
            parts = [int(x) for x in s.replace(",", " ").split()]
        except ValueError:
            return "must be 4 integers"
        if len(parts) != 4:
            return "must be exactly 4 integers (D H M S)"
        if parts == [0, 0, 0, 0]:
            return "must not be all zero"
        return None

    def source_start(v, val):
        return fourvar(v, val) or hms_match(v, val)

    def secflds(v, val):
        from misc import secflds_list
        if isinstance(val, str) and val.strip() == "":
            return None
        try:
            secflds_list(val)
        except ValueError as e:
            return str(e)
        return None

    def other_input(v, val):                        # keys must be &tgcm_input members
        from misc import as_list, unknown_other_input_keys, other_input_error
        lines = val if isinstance(val, (list, tuple)) else as_list(val)
        unknown = unknown_other_input_keys(lines)
        return other_input_error(unknown) if unknown else None

    reg = {
        "inp.SECFLDS": secflds,
        "inp.other_input": other_input,
        "inp.start_time": iso,
        "inp.stop_time": stop_after,
        "inp.SOURCE_START": source_start,
        "inp.segment": seg_required,
        "inp.secondary_start_time": sec_window("secondary_start_time"),
        "inp.secondary_stop_time": sec_window("secondary_stop_time"),
    }
    for fld in ("PRISTART", "PRISTOP", "PRIHIST", "SECHIST", "SECSTART", "SECSTOP"):
        reg["inp." + fld] = fourvar

    def forcing(fld):                               # required by input.F without GPI
        def check(v, val):
            from namelist_solver import forcing_missing
            inp = {p[4:]: x for p, x in v.items() if p.startswith("inp.")}
            return dict(forcing_missing(inp)).get(fld)
        return check
    for fld in ("POWER", "CTPOTEN", "F107", "F107A", "BXIMF", "BYIMF", "BZIMF", "SWDEN", "SWVEL"):
        reg["inp." + fld] = forcing(fld)

    def account(machine):
        def check(v, val):
            if _is_set(val):
                return None
            return f"required on {machine}: the PBS account (#PBS -A) the jobs are charged to"
        return check
    from output_solver import CONFIG_DIR
    from jobgen import JobGen
    for name, m in JobGen(CONFIG_DIR, machine="derecho").machines.items():
        if m.get("scheduler", "none") != "none":
            reg[f"job.{name}.project_code"] = account(name)
    return reg


def _warn_registry():
    """Return {path: fn(values) -> warning or None}, the warnings that depend on other fields."""
    from namelist_solver import inp_pri_date, inp_mxhist
    from misc import segment_time

    _gpy_missing = {}

    def _gcmprocpy_missing():
        if "v" not in _gpy_missing:
            import importlib.util
            _gpy_missing["v"] = importlib.util.find_spec("gcmprocpy") is None
        return _gpy_missing["v"]

    def seg(v):
        s = str(v["inp.segment"]).strip()
        if s in ("", "[none]", "[None]", "[null]", "none", "None", "[]"):
            return None
        try:
            out = [int(x) for x in s.split()]
        except ValueError:
            return None
        return out if any(out) else None

    def seg_label(label):
        def fn(v):
            s = seg(v)
            if s is None:
                return None
            rt = segment_time(v["inp.start_time"], v["inp.stop_time"], s)
            return (f"Segmentation is set to {s}.\n{label} is set for one segment"
                    f"\neg.{rt[0][0]} to {rt[0][1]}")
        return fn

    def mxhist(hist_path):
        def fn(v):
            xh = [int(x) for x in str(v[hist_path]).split()]
            return inp_mxhist(v["inp.start_time"], v["inp.stop_time"], xh, None, seg(v))[1]
        return fn

    def year_boundary(v):
        if v.get("__benchmark"):
            return None
        pri = inp_pri_date(v["inp.start_time"], v["inp.stop_time"])
        year, pristop_day = int(pri[0]), int(pri[3][0])
        leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
        mxday = 367 if leap else 366
        if pristop_day <= mxday:
            return None
        return (f"This run crosses the {year}->{year + 1} year boundary (it would run to day "
                f"{pristop_day}, but the model caps a single run at day {mxday}). Segmentation is "
                f"required: the run is split automatically at Jan 1. Enter a segment length, "
                f"e.g. '5 0 0 0'.")

    def if_seg(text):
        return lambda v: (text if seg(v) is not None else None)

    def gpi_gen(v):
        if v.get("__benchmark"):
            return None
        if str(v["inp.GPI_NCFILE"]).strip().lower() == "gen":
            if _gcmprocpy_missing():
                return ("'gen' → a GPI file would be generated for this run's dates on accept, but "
                        "gcmprocpy is not installed in this environment, so generation would fail. "
                        "Install it (pip install gcmprocpy) or enter a GPI file path / 'none'.")
            return ("'gen' → a GPI file will be generated for this run's dates on accept, with "
                    "gcmprocpy (27-day trailing F10.7 average).")
        return ("Enter 'gen' to generate a fresh GPI file for this run's dates with gcmprocpy "
                "(27-day trailing F10.7 average) instead of this one.")

    def imf_gen(v):
        if v.get("__benchmark"):
            return None
        if str(v["inp.IMF_NCFILE"]).strip().lower() == "gen" and _gcmprocpy_missing():
            return ("'gen' → an IMF file would be generated on accept, but gcmprocpy is not "
                    "installed here, so generation would fail. Install it (pip install gcmprocpy) "
                    "or enter an IMF file path / 'none' to skip.")
        return ("Enter 'gen' to generate an IMF file for this run's dates with gcmprocpy "
                "(OMNI solar wind).")

    def under_gpi(text):
        return lambda v: (text if _is_set(v["inp.GPI_NCFILE"]) else None)

    def exe_warn(path):
        def fn(v):
            if v.get("__compile") or v.get("__onlycompile"):
                return None
            exe = str(v.get(path, "") or "").strip()
            if exe in ("", "None", "none", "[None]"):
                return None
            if os.path.isfile(exe):
                return None
            execdir = str(v.get("model.data.execdir", "") or "")
            if execdir and os.path.isfile(os.path.join(execdir, os.path.basename(exe))):
                return None
            return (f"{exe} not found — the model must be compiled before it runs. Re-run with "
                    f"--compile/-c (build then run) or --onlycompile/-oc (build only).")
        return fn

    derives = _derive_registry()

    def derived(path, then=None):
        def fn(v):
            notes = []
            try:
                want = derives[path][0](v)
                cur = v.get(path)
                if str(cur).strip() not in ("", "None") and str(cur).strip() != str(want).strip():
                    notes.append("not used: a segmented run derives it per segment"
                                 if _form_segmented(v) else f"differs from the derived {want}")
            except Exception:
                pass
            extra = then(v) if then else None
            return "\n".join(notes + ([extra] if extra else [])) or None
        return fn

    return {
        "inp.segment":     year_boundary,
        "inp.PRIHIST":     seg_label("PRIHIST"),
        "inp.SECHIST":     seg_label("SECHIST"),
        "inp.MXHIST_PRIM": mxhist("inp.PRIHIST"),
        "inp.MXHIST_SECH": mxhist("inp.SECHIST"),
        **{f"inp.{k}": derived(f"inp.{k}") for k in ("START_YEAR", "START_DAY", "PRISTART", "PRISTOP",
                                                    "SECSTART", "SECSTOP")},
        "inp.OUTPUT":      derived("inp.OUTPUT", if_seg("Primary Output can be ignored. Will be set on segmentation")),
        "inp.SECOUT":      derived("inp.SECOUT", if_seg("Secondary Output can be ignored. Will be set on segmentation")),
        "inp.GPI_NCFILE":  gpi_gen,
        "inp.IMF_NCFILE":  imf_gen,
        "inp.F107":        under_gpi("F10.7 can be read by GPI File and can be skipped."),
        "inp.F107A":       under_gpi("81-Day Average of F10.7 can be read by GPI File and can be skipped."),
        "model.data.modelexe":         exe_warn("model.data.modelexe"),
        "model.data.coupled_modelexe": exe_warn("model.data.coupled_modelexe"),
    }


def _inp(v):
    """Return the inp.* values as an inp dict."""
    return {p[4:]: x for p, x in v.items() if p.startswith("inp.")}


def _is_set(x):
    """True unless x is blank, null, false or 0."""
    s = str(x).strip().lower()
    return s not in ("", "none", "null", "[null]", "[none]", "false", "0", "[]")


def _cond_registry():
    """Return cond(field, values) -> True to show; hides fields the run would ignore."""
    def cond(f, v):
        p = f.path
        hpc = str(v.get("simulation.hpc_system", "")).strip().lower()
        bench = bool(v.get("__benchmark"))
        coupling = bool(v.get("__coupling"))

        if f.name in v.get("__engage_skip", ()):
            return False
        if p == "model.data.input_file" and (bench or v.get("__engage")):
            return False

        if v.get("__onlycompile"):
            if p.startswith("inp.") or p.startswith("job."):
                return False
            if p in ("model.data.input_file", "model.data.log_file", "simulation.job_name",
                     "model.data.modeldir", "model.data.parentdir", "model.data.tgcmdata",
                     "model.specification.segmentation"):
                return False
            if coupling and p == "model.data.modelexe":
                return False

        # a custom input_file supplies the whole inp namelist
        if p.startswith("inp.") and _is_set(v.get("model.data.input_file")):
            return False
        if p == "model.data.coupled_modelexe" and not coupling:
            return False
        if p.startswith("job."):
            if hpc == "linux":
                return False
            if p.split(".")[1] != hpc:
                return False
            # engage_run derives the coupled PE layout
            if v.get("__engage") and (".resource." in p or p.endswith(".nprocs")):
                return False
        if v.get("__engage"):
            if p in ("model.data.execdir", "model.data.workdir", "model.data.histdir",
                     "model.specification.nres_grid"):
                return False
        if bench and p in ("inp.start_time", "inp.stop_time", "inp.secondary_start_time",
                           "inp.secondary_stop_time", "inp.segment", "inp.solar_flux_level"):
            return False
        pot = str(v.get("inp.POTENTIAL_MODEL", "HEELIS")).strip().upper()
        gpi = _is_set(v.get("inp.GPI_NCFILE"))
        imf = _is_set(v.get("inp.IMF_NCFILE"))
        kp = _is_set(v.get("inp.KP"))
        if p in ("inp.IMF_NCFILE", "inp.BXIMF", "inp.BYIMF", "inp.BZIMF", "inp.SWDEN",
                 "inp.SWVEL") and pot != "WEIMER":
            return False
        if p in ("inp.BXIMF", "inp.BYIMF", "inp.BZIMF", "inp.SWDEN", "inp.SWVEL") and imf:
            return False
        if p in ("inp.KP", "inp.POWER", "inp.CTPOTEN", "inp.F107", "inp.F107A") and gpi:
            return False
        if p in ("inp.POWER", "inp.CTPOTEN") and kp:
            return False
        if p in ("inp.CTPOTEN", "inp.KP") and pot == "WEIMER":
            return False
        if f.name in _GSWM_NM_FIELDS:
            try:
                if float(v.get("model.specification.horires")) == 5.0:
                    return False
            except (TypeError, ValueError):
                pass
        return True
    return cond


def build_state(args, mode="EXPERT"):
    """Build the FormState from options_description.json and the run context."""
    here = os.path.dirname(os.path.abspath(__file__))
    from output_solver import apply_system_config
    with open(os.path.join(here, "options_description.json"), encoding="utf-8") as f:
        od = apply_system_config(json.load(f), socket.gethostname())
    fields = load_tiegcmrun_fields(od)
    bench = getattr(args, "benchmark", None)
    eng = getattr(args, "engage", None)
    eng = eng if isinstance(eng, dict) else None
    bench_inp = {}
    if bench:
        try:
            with open(os.path.join(here, "benchmarks.json"), encoding="utf-8") as f:
                bench_inp = (json.load(f).get(bench, {}) or {}).get("inp", {}) or {}
        except (OSError, ValueError):
            bench_inp = {}
    ctx = {"__engage": eng is not None,
           "__engage_cfg": eng,
           "__engage_skip": set((eng or {}).get("skip", ())),
           "__benchmark": bench is not None,
           "__benchmark_name": bench,
           "__bench_inp": bench_inp,
           "__coupling": getattr(args, "coupling", False),
           "__onlycompile": getattr(args, "onlycompile", False),
           "__compile": getattr(args, "compile", False)}
    if eng is not None:
        from engage_solver import COUPLED_DEFAULT_KEYS
        for f in fields:
            if f.section == "inp" and f.name in COUPLED_DEFAULT_KEYS:
                f.level = "EXPERT"
    state = FormState(fields, _derive_registry(), _validate_registry(), mode=mode,
                      context=ctx, cond=_cond_registry(), warn=_warn_registry())
    _seed_benchmark_engage(state)
    return state


def _seed_benchmark_engage(state):
    """Pin the engage run window/grid and the benchmarks.json inp overrides."""
    from misc import resolve_benchmark_file
    ctx = state.context
    tg = tiegcm_env("TIEGCMDATA") or ""
    eng = ctx.get("__engage_cfg")
    if eng:
        for key, path in (("start_time", "inp.start_time"), ("stop_time", "inp.stop_time"),
                          ("segment", "inp.segment"), ("horires", "model.specification.horires")):
            if eng.get(key) is not None:
                val = eng[key]
                # [0,2,0,0] -> '0 2 0 0'
                if key == "segment" and isinstance(val, list):
                    val = " ".join(map(str, val))
                state.set(path, val)
    bench_inp = ctx.get("__bench_inp") or {}
    if bench_inp:
        jn = state.values["simulation.job_name"]
        run_name = f'{jn}_{state.values["model.specification.horires"]}x{state.values["model.specification.vertres"]}'
        from misc import history_dir
        histdir = history_dir({"parentdir": state.values.get("model.data.parentdir"),
                               "workdir": state.values["model.data.workdir"],
                               "histdir": state.values["model.data.histdir"]})
        for fld, raw in bench_inp.items():
            path = "inp." + fld
            if path not in state.by_path:
                continue
            if raw is None:
                # null keeps the static default; a derived field is pinned so its derive is unused
                if path in state.derive:
                    dflt = state.by_path[path].default
                    state.set(path, "" if dflt is None else dflt)
                continue
            if fld == "SOURCE_START":
                continue        # derived from SOURCE's histories
            if fld in ("SOURCE", "GPI_NCFILE", "IMF_NCFILE"):
                val = resolve_benchmark_file(fld, ctx.get("__benchmark_name"), bench_inp, tg)
            elif fld in ("OUTPUT", "SECOUT"):
                val = str(raw).replace("+histdir+", histdir).replace("+run_name+", run_name)
            elif fld == "other_input":
                val = [it.replace("+tiegcmdata+", tg) if isinstance(it, str) else it for it in raw]
            elif (fld in ("SOURCE_START", "segment", "PRISTART", "PRISTOP", "PRIHIST",
                          "SECSTART", "SECSTOP", "SECHIST") and isinstance(raw, list)):
                val = " ".join(map(str, raw))
            else:
                val = raw
            state.set(path, val)


def run_tui_form(args):
    """Run the all-fields form and return the options dict; _FormUnavailable without a TTY."""
    import sys
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise _FormUnavailable("no TTY")
    try:
        from prompt_toolkit import Application
        from prompt_toolkit.key_binding import KeyBindings
        from prompt_toolkit.layout import Layout
        from prompt_toolkit.layout.containers import HSplit, VSplit, Window
        try:
            from prompt_toolkit.layout import ScrollablePane
        except ImportError:
            ScrollablePane = None                                        # old prompt_toolkit: no scrolling
        from prompt_toolkit.layout.controls import FormattedTextControl
        from prompt_toolkit.layout.dimension import D
        from prompt_toolkit.widgets import TextArea, Frame
        from prompt_toolkit.styles import Style
    except ImportError as e:
        raise _FormUnavailable(str(e))

    state = build_state(args, mode=(getattr(args, "mode", None) or "EXPERT"))
    vis = state.visible()
    prog = {"on": False}
    areas = {}
    for f in vis:
        # min=10 so a long default cannot push the form past the terminal width
        areas[f.path] = TextArea(text=str(state.values[f.path]), multiline=False, height=1,
                                 width=D(weight=1, min=10))
    err = FormattedTextControl(text="")

    def refresh():
        prog["on"] = True
        for p, ta in areas.items():
            cur = str(state.values[p])
            if ta.text != cur:
                ta.text = cur
        prog["on"] = False
        es = state.errors()
        err.text = ([("class:ok", "  ✓ valid — Ctrl-S to accept, Ctrl-C to cancel\n")] if not es
                    else [("class:err", f"  ✗ {p.split('.')[-1]}: {m}\n") for p, m in list(es.items())[:6]])

    def on_change(path):
        def _h(buf):
            if not prog["on"]:
                state.set(path, buf.text, normalize=False)   # per keystroke: no clean-up mid-typing
                refresh()
        return _h

    rows = []
    last_section = None
    for f in vis:
        areas[f.path].buffer.on_text_changed += on_change(f.path)
        if f.section != last_section:
            rows.append(Window(FormattedTextControl([("class:section", f" {f.section}")]), height=1))
            last_section = f.section
        tag = "*" if f.path in state.derive else " "
        rows.append(VSplit([
            Window(FormattedTextControl(f" {f.name:<22}"), width=24),
            areas[f.path],
            Window(FormattedTextControl(f"{tag}[{f.level[0]}]"), width=5),
        ], height=1))

    field_pane = HSplit(rows)
    title = Window(FormattedTextControl(
        [("class:title", f" tiegcmrun [{state.mode}]  Tab/arrows: move - type: edit - Ctrl-S: accept - Ctrl-C: cancel ")]),
        height=1, style="class:title")
    body = HSplit([
        title,
        ScrollablePane(field_pane) if ScrollablePane else field_pane,
        Window(content=err, height=4),
    ])
    root = body

    kb = KeyBindings()
    result = {"opts": None}

    @kb.add("c-c")
    def _(e):
        e.app.exit()

    @kb.add("c-s")
    def _(e):
        if not state.errors():
            result["opts"] = state.options_dict()
            e.app.exit()

    @kb.add("down")
    @kb.add("tab")
    def _(e):
        e.app.layout.focus_next()

    @kb.add("up")
    @kb.add("s-tab")
    def _(e):
        e.app.layout.focus_previous()

    style = Style.from_dict({"ok": "ansigreen", "err": "ansired bold",
                             "section": "bold ansicyan", "title": "reverse"})
    # no mouse_support or focused_element: both cause terminal-specific render errors
    app = Application(layout=Layout(root), key_bindings=kb, style=style, full_screen=True)
    refresh()
    app.run()
    if result["opts"] is None:
        raise _FormCancelled("cancelled")
    return result["opts"]


def run_tui_wizard(args):
    """Ask the options one at a time with arrow navigation; return the options dict.

    Going back and editing re-derives the later defaults that were not typed.
    """
    import sys
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise _FormUnavailable("no TTY")
    try:
        from prompt_toolkit import Application
        from prompt_toolkit.application import get_app
        from prompt_toolkit.key_binding import KeyBindings
        from prompt_toolkit.filters import Condition
        from prompt_toolkit.layout import Layout
        from prompt_toolkit.layout.containers import ConditionalContainer, HSplit, VSplit, Window
        from prompt_toolkit.layout.controls import FormattedTextControl
        from prompt_toolkit.layout.dimension import D
        from prompt_toolkit.widgets import TextArea
        from prompt_toolkit.styles import Style
        try:
            from prompt_toolkit.layout import ScrollablePane
        except ImportError:
            ScrollablePane = None
    except ImportError as e:
        raise _FormUnavailable(str(e))

    mode = getattr(args, "mode", None)
    if mode is None:
        mode = "BENCH" if getattr(args, "benchmark", None) is not None else "BASIC"
    state = build_state(args, mode=mode)
    vis0 = state.visible()
    if not vis0:
        raise _FormUnavailable("no visible fields")
    cur = {"path": vis0[0].path}
    blocked = {"on": False}                       # a forward move was refused
    result = {"opts": None}
    gen_failed = {}                               # path -> why a GPI/IMF 'gen' failed
    accept_failed = {}                            # path -> why the accepted options were refused

    def vis_paths():
        return [f.path for f in state.visible()]

    input_area = TextArea(text=str(state.values[cur["path"]]), multiline=False, height=1)
    input_area.buffer.cursor_position = len(input_area.text)

    def load(path):
        cur["path"] = path
        input_area.text = str(state.values[path])
        input_area.buffer.cursor_position = len(input_area.text)

    def commit():
        p = cur["path"]
        typed = input_area.text
        if typed != str(state.values[p]):       # only an actual edit pins the value
            state.set(p, typed)
            gen_failed.pop(p, None)
            accept_failed.clear()

    def picks(p):
        # valids, else suggestions ('choices': free text still accepted)
        f = state.by_path[p]
        return [str(x) for x in (f.valids or (f.meta or {}).get("choices") or [])]

    def cycle_valid(delta):
        p = cur["path"]
        valids = picks(p)
        if not valids:
            return
        cur_text = input_area.text.strip()
        try:
            i = valids.index(cur_text)
        except ValueError:
            i = -1 if delta > 0 else 0           # first press lands on valids[0] / valids[-1]
        i = (i + delta) % len(valids)
        input_area.text = valids[i]
        input_area.buffer.cursor_position = len(input_area.text)
        state.set(p, valids[i])
        accept_failed.clear()

    def finish():
        commit()
        errs = state.errors()
        if not errs:
            result["opts"] = state.options_dict()
            app.exit()
        else:
            for p in vis_paths():
                if p in errs:
                    load(p)
                    return

    def move(delta):
        commit()
        if delta > 0 and cur["path"] in state.errors():
            blocked["on"] = True
            return
        blocked["on"] = False
        vp = vis_paths()
        p = cur["path"]
        if p not in vp:                          # the current field hid itself
            load(vp[0])
            return
        j = vp.index(p) + delta
        if j < 0:
            load(vp[0])
        elif j >= len(vp):
            finish()
        else:
            load(vp[j])

    def tag(p):
        if p in state.overridden:
            return "your value"
        return "derived" if p in state.derive else "default"

    def header():
        vp = vis_paths()
        p = cur["path"]
        i = (vp.index(p) + 1) if p in vp else 0
        f = state.by_path[p]
        return [("class:title", f" tiegcmrun wizard — {state.mode}    step {i}/{len(vp)}    "
                                f"[{f.section}] ")]

    def _shown(raw):
        return "—" if (raw is None or str(raw).strip() in ("", "[None]", "[none]", "None",
                                                           "[null]")) else str(raw)

    def overview():
        vp = vis_paths()
        p = cur["path"]
        ci = vp.index(p) if p in vp else 0
        lo, hi = max(0, ci - 9), min(len(vp), ci + 3)
        out = []
        if lo > 0:
            out.append(("class:more", f"     … {lo} earlier\n"))
        for i in range(lo, hi):
            f = state.by_path[vp[i]]
            val = _shown(state.values[vp[i]])
            if len(val) > 48:
                val = val[:47] + "…"
            if i == ci:
                out.append(("class:cur", f"  ▶ {f.name:<22} {val}"))
                out.append(("class:cur", "\n"))
            else:
                mark = "✓" if vp[i] in state.overridden else " "
                out.append(("class:donemark", f"  {mark} "))
                out.append(("class:donename", f"{f.name:<22} "))
                out.append(("class:doneval", f"{val}\n"))
        if hi < len(vp):
            out.append(("class:more", f"     … {len(vp) - hi} more\n"))
        return out

    def body():
        p = cur["path"]
        f = state.by_path[p]
        out = [("class:name", f"  {f.name}"), ("class:leveltag", f"   [{f.level}]\n")]
        for ln in str(f.prompt).split("\n"):
            out.append(("class:prompttext", f"  {ln}\n"))
        if f.description:
            out.append(("class:help", f"  {f.description}\n"))
        shown = picks(p)
        if shown:
            curv = str(state.values[p]).strip()
            chips = "  ".join((f"[{v}]" if str(v) == curv else f" {v} ") for v in shown)
            how = "(←/→ to choose)" if f.valids else "(←/→ to choose, or type another)"
            out.append(("class:valids", f"  {chips}   {how}\n"))
        out.append(("class:dim", f"  default ({tag(p)}):  "))
        out.append(("class:default", f"{_shown(state.values[p])}\n"))
        if f.warning:
            out.append(("class:warn", f"  ⚠ {f.warning}\n"))
        dw = state.warnings_for(p, live=(input_area.text if p in _LIVE_WARN else None))
        if dw and dw != f.warning:
            for ln in str(dw).split("\n"):
                out.append(("class:warn", f"  ⚠ {ln}\n"))
        if p in state.derive_errors:
            out.append(("class:warn", f"  (could not derive — check upstream: "
                                      f"{state.derive_errors[p]})\n"))
        return out

    def errbox():
        # own window above the input row, so a long help text cannot clip the error
        p = cur["path"]
        out = []
        for msgs in (gen_failed, accept_failed):
            if p in msgs:
                for ln in str(msgs[p]).split("\n"):
                    out.append(("class:err", f"  ✗ {ln}\n"))
        err = state.errors().get(p)
        if err:
            out.append(("class:err", f"  ✗ {err}\n"))
            if blocked["on"]:
                out.append(("class:err", "  ↳ fix this to continue (or ↑ to go back)\n"))
        if out:
            out[-1] = (out[-1][0], out[-1][1].rstrip("\n"))
        return out

    def footer():
        # status first: a narrow terminal truncates the end of the line
        n = len(state.errors())
        try:
            cols = get_app().output.get_size().columns
        except Exception:
            cols = 120
        wide = cols >= 100
        status = (("class:err", f"  ({n} field(s) need fixing)") if n
                  else ("class:ok", "  (all valid — Ctrl-S to finish)" if wide
                        else "  (all valid)"))
        legend = ("   ↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel"
                  if wide else "  ↑ back  ↓/Enter next  ^S finish  ^R reset  ^C cancel")
        return [status, ("class:foot", legend)]

    left = HSplit([
        Window(FormattedTextControl(body), dont_extend_height=True, wrap_lines=True),
        ConditionalContainer(
            Window(FormattedTextControl(errbox), dont_extend_height=True, wrap_lines=True),
            filter=Condition(lambda: bool(errbox()))),
        VSplit([Window(FormattedTextControl([("class:caret", "  › ")]), width=4), input_area]),
        Window(),
    ], width=D(weight=3))
    overview_win = Window(FormattedTextControl(overview), dont_extend_height=True)
    right = (ScrollablePane(overview_win, width=D(weight=2)) if ScrollablePane
             else HSplit([overview_win], width=D(weight=2)))
    root = HSplit([
        Window(FormattedTextControl(header), height=1, style="class:title"),
        VSplit([
            left,
            Window(width=1, char="│", style="class:rule"),
            right,
        ]),
        Window(FormattedTextControl(footer), height=1, style="class:title"),
    ])

    kb = KeyBindings()

    @kb.add("up")
    def _(e):
        move(-1)

    @kb.add("down")
    @kb.add("enter")
    def _(e):
        move(1)

    # ←/→ cycle the valids; free-text fields keep normal cursor movement
    _on_valids = Condition(lambda: bool(picks(cur["path"])))

    @kb.add("left", filter=_on_valids)
    def _(e):
        cycle_valid(-1)

    @kb.add("right", filter=_on_valids)
    def _(e):
        cycle_valid(1)

    @kb.add("c-s")
    def _(e):
        finish()

    @kb.add("c-r")
    def _(e):
        p = cur["path"]
        state.reset(p)
        gen_failed.pop(p, None)
        accept_failed.clear()
        blocked["on"] = False
        load(p)

    @kb.add("c-c")
    def _(e):
        e.app.exit()

    style = Style.from_dict({
        "title":     "reverse",
        "cur":       "reverse bold",
        "donemark":  "ansigreen",
        "donename":  "ansibrightblack",
        "doneval":   "ansicyan",
        "more":      "ansibrightblack italic",
        "rule":      "ansibrightblack",
        "name":      "bold ansicyan",
        "leveltag":  "ansibrightblack",
        "prompttext": "",
        "help":      "ansibrightblack",
        "valids":    "ansimagenta",
        "dim":       "ansibrightblack",
        "default":   "ansigreen bold",
        "warn":      "ansiyellow",
        "err":       "ansired bold",
        "caret":     "ansigreen bold",
        "foot":      "ansibrightblack",
        "ok":        "ansigreen",
    })
    # no mouse_support or focused_element: both cause terminal-specific render errors
    app = Application(layout=Layout(root), key_bindings=kb, style=style, full_screen=True)
    while True:
        app.run()
        if result["opts"] is None:
            raise _FormCancelled("cancelled")
        # A failed 'gen' or a refused check re-opens the wizard on the field it names.
        drop_segmented_derived(state)
        failed = generate_pending(state)
        if not failed:
            failed = accept_problems(state, args)
            if not failed:
                note_coupled_unset(state, args)
                return state.options_dict()
            accept_failed.update(failed)
        else:
            gen_failed.update(failed)
        result["opts"] = None
        blocked["on"] = False
        vp = vis_paths()
        load(next((p for p in vp if p in failed), vp[0]))


_PROBLEM_KEY = re.compile(r"\b((?:inp|job|model\.data|model\.specification|simulation)(?:\.\w+)+)")


def _problem_field(problem, paths):
    """Return the field a validate_options problem names, else the last field."""
    m = _PROBLEM_KEY.search(problem)
    if m:
        key = m.group(1)
        for p in paths:
            parts = p.split(".")
            norm = ".".join(parts[:1] + parts[2:]) if parts[0] == "job" else p
            if norm == key or norm.startswith(key + ".") or key.startswith(norm + "."):
                return p
    return paths[-1]


def accept_problems(state, args):
    """Run tiegcmrun's pre-write validation on the form's options; return {path: problems}."""
    import contextlib
    import io
    from output_solver import apply_system_config
    from replay_solver import validate_options
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "options_description.json"), encoding="utf-8") as f:
        od = apply_system_config(json.load(f), socket.gethostname())
    probe = finalize_options(copy.deepcopy(state.options_dict()), args)
    with contextlib.redirect_stdout(io.StringIO()):     # NOTEs print again on the real run
        problems = validate_options(probe, od)
    out = {}
    paths = [f.path for f in state.visible()]
    for problem in problems:
        p = _problem_field(problem, paths)
        out[p] = f"{out[p]}\n{problem}" if p in out else problem
    return out


_DERIVED_PATHS = tuple(f"inp.{k}" for k in ("START_YEAR", "START_DAY", "PRISTART", "PRISTOP",
                                             "SECSTART", "SECSTOP", "OUTPUT", "SECOUT"))


def _form_segmented(v):
    """True when the run is segmented (an engage run always is)."""
    if v.get("__engage"):
        return True
    s = str(v.get("inp.segment")).strip()
    try:
        return any(int(x) for x in s.split())
    except ValueError:
        return False


def drop_segmented_derived(state):
    """Reset edited run-window keys of a segmented run, which derives them per segment."""
    if not _form_segmented(state._all()):
        return
    for path in _DERIVED_PATHS:
        if path not in state.overridden or path not in state.derive:
            continue
        cur = state.values.get(path)
        state.reset(path)
        if str(cur).strip() not in ("", "None") and str(cur).strip() != str(state.values.get(path)).strip():
            print(f"\033[33mNOTE: {path[4:]} {cur} not used: a segmented run derives it per segment.\033[0m")


def note_coupled_unset(state, args):
    """Record the coupled default keys the user cleared, so engage leaves them unset."""
    eng = getattr(args, "engage", None)
    if not isinstance(eng, dict):
        return
    from engage_solver import COUPLED_DEFAULT_KEYS
    eng["coupled_unset"] = [k for k in COUPLED_DEFAULT_KEYS if f"inp.{k}" in state.overridden
                            and not _is_set(state.values.get(f"inp.{k}"))]


def _gen_fields():
    import misc
    return (("inp.GPI_NCFILE", "GPI", misc.generate_gpi_file),
            ("inp.IMF_NCFILE", "IMF", misc.generate_imf_file))


def generate_pending(state):
    """Generate each visible GPI/IMF field set to 'gen'; return {path: message} for failures.

    A failed field is blanked so the forcing fields it replaced become required.
    """
    import misc
    tg = state.values.get("model.data.tgcmdata")
    if _is_set(tg):
        os.environ["TIEGCMDATA"] = tg
    start, stop = state.values.get("inp.start_time"), state.values.get("inp.stop_time")
    workdir = state.values.get("model.data.workdir")
    visible = {f.path for f in state.visible()}
    failed = {}
    for path, kind, gen in _gen_fields():
        val = state.values.get(path)
        if path not in visible or not (isinstance(val, str) and val.strip().lower() == "gen"):
            continue
        print(f"==> {path[4:]}: generating a {kind} file for {start}..{stop} into {workdir} ...",
              flush=True)
        try:
            out = gen(start, stop, workdir)
        except Exception as e:
            print(f"{path[4:]}: generation failed: {e}", flush=True)
            state.set(path, "")
            failed[path] = (f"{kind} generation failed: {e}\n"
                            f"Enter a {kind} file path (or 'gen' to retry), or leave it blank and "
                            f"give the forcing values instead.")
            continue
        print(f"==> {path[4:]}: generated {out}", flush=True)
        if not misc.file_covers(out, start, stop):
            print(f"WARNING {path[4:]}: generated file {out} does not span the full run "
                  f"({start}..{stop}); the indices source may lag real time.", flush=True)
        state.set(path, out)
    return failed


def apply_accept_actions(options, args):
    """Apply the linear prompts' inline side effects to the accepted options; return them.

    Sets TIEGCMDATA, generates any remaining GPI/IMF 'gen' (raises on failure) and resolves exes.
    """
    import misc
    data = options.get("model", {}).get("data", {})
    inp = options.get("inp", {})

    tg = data.get("tgcmdata")
    if tg:
        os.environ["TIEGCMDATA"] = tg

    misc.generate_pending(inp, data.get("workdir"))

    nocompile = not getattr(args, "compile", False) and not getattr(args, "onlycompile", False)

    def _exe(key):
        exe = data.get(key)
        if not exe or os.path.isfile(exe):
            return
        exe = misc.in_execdir(data.get("execdir", ""), exe)
        data[key] = exe
        if nocompile and not os.path.isfile(exe):
            print(f"Warning: {exe} not found — the model must be compiled (--compile / --onlycompile).")
    _exe("modelexe")
    if getattr(args, "coupling", False):
        _exe("coupled_modelexe")
    return finalize_options(options, args)


def finalize_options(options, args):
    """Reshape the form's options as the linear prompts build them; return them (no side effects)."""
    data = options.get("model", {}).get("data", {})
    onlycompile = getattr(args, "onlycompile", False)
    if _is_set(data.get("input_file")) or onlycompile:
        options.pop("inp", None)
    if data.get("modeldir"):
        data["utildir"] = os.path.join(data["modeldir"], "scripts")
    if not getattr(args, "coupling", False):
        data.pop("coupled_modelexe", None)
    # job = _common + the active machine's block
    job = options.get("job")
    if isinstance(job, dict):
        hpc = str(options.get("simulation", {}).get("hpc_system", "")).strip()
        if hpc == "linux" or onlycompile:
            options.pop("job", None)
        else:
            flat = {}
            flat.update(job.get("_common", {}) or {})
            flat.update(job.get(hpc, {}) or {})
            options["job"] = flat
    return options
