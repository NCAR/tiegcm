"""Namelist value helpers for tiegcmrun: run dates, history cadences and files, forcing and
run-window checks."""


import re
from datetime import datetime, timedelta
from math import ceil, gcd

RUN_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S"


def parse_run_datetime(value, field="date"):
    """Parse YYYY-MM-DDThh:mm:ss strictly; ValueError naming field otherwise."""
    if isinstance(value, datetime):
        return value
    try:
        return datetime.strptime(str(value).strip(), RUN_DATETIME_FORMAT)
    except (TypeError, ValueError):
        raise ValueError(f"{field} must be given as yyyy-mm-ddThh:mm:ss "
                         f"(e.g. 2013-03-17T00:00:00); got {value!r}") from None


def _mtime_seconds(mtime):
    d, h, m, s = [int(x) for x in mtime]
    return d * 86400 + h * 3600 + m * 60 + s


def _seconds_mtime(seconds):
    seconds = int(seconds)
    return [seconds // 86400, (seconds % 86400) // 3600, (seconds % 3600) // 60, seconds % 60]


def parse_segment(value):
    """A segment length 'D H M S' as [d,h,m,s], or None for unset or all zeros."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if len(value) == 0 or all(x is None for x in value):
            return None
        parts = list(value)
    else:
        s = str(value).strip()
        if s.lower() in ("", "none", "null", "[none]", "[null]", "[]"):
            return None
        parts = s.replace(",", " ").split()
    try:
        seg = [int(x) for x in parts]
    except (TypeError, ValueError):
        raise ValueError(f"segment must be 4 integers 'D H M S' (e.g. '5 0 0 0'); got {value!r}") from None
    if len(seg) != 4:
        raise ValueError(f"segment must be 4 integers 'D H M S' (e.g. '5 0 0 0'); got {value!r}")
    return seg if any(seg) else None


def fit_cadence(cadence, span_sec, step=60):
    """The history cadence [d,h,m,s] nearest to cadence that divides span_sec and is a multiple
    of step (input.F requires both); ValueError when span_sec is not a multiple of step."""
    step = int(step)
    if step <= 0:
        raise ValueError(f"fit_cadence: STEP must be a positive number of seconds; got {step}")
    c = _mtime_seconds(cadence)
    span_sec = int(span_sec)
    if c <= 0 or span_sec <= 0:
        return cadence
    if span_sec % step != 0:
        raise ValueError(
            f"The history window is {span_sec} s, which is not a multiple of STEP {step} s: TIE-GCM "
            f"writes histories only on the STEP clock, so no history cadence fits it. Start and stop "
            f"the run on a multiple of STEP.")
    if c % step != 0:
        c = -(-c // step) * step
        cadence = _seconds_mtime(c)
    if span_sec % c == 0:
        return cadence
    if span_sec < c:
        return _seconds_mtime(span_sec)
    return _seconds_mtime(gcd(span_sec, c))


def cadence_step(engage, *steps):
    """The STEP (s) a default cadence must be a multiple of: engage's, else the first set step."""
    if engage is not None and engage.get("STEP") is not None:
        return int(engage["STEP"])
    for step in steps:
        if step is not None and str(step).strip().lower() not in ("", "none", "null"):
            try:
                return int(float(step))
            except (TypeError, ValueError):
                continue
    return 60


def segment_mxhist(n_hist, mxhist):
    """MXHIST_PRIM for a segment of n_hist histories: the largest divisor of n_hist <= mxhist,
    so the segment's last history ends up alone in the file the next segment reads as SOURCE."""
    mx = int(mxhist)
    n_hist = int(n_hist)
    if n_hist < 1 or mx < 1:
        return mx
    mx = min(mx, n_hist)
    while n_hist % mx:
        mx -= 1
    return mx


def nc_wrhist_mtime(mtime, start_year):
    """The mtime nc_wrhist stores for a history at model time mtime: the day wraps after the end
    of start_year. A continuing segment must use this as SOURCE_START."""
    mt = [int(x) for x in mtime]
    leap = (start_year % 4 == 0 and start_year % 100 != 0) or start_year % 400 == 0
    mfinal = 366 if leap else 365
    if mt[0] > mfinal:
        mt[0] -= mfinal
    return mt


def _segment_cadence(segment):
    """One unit of the segment's finest non-zero field, which always divides the segment."""
    if segment[3] != 0:
        return [0, 0, 0, 1]
    elif segment[2] != 0:
        return [0, 0, 1, 0]
    elif segment[1] != 0:
        return [0, 1, 0, 0]
    else:
        return [1, 0, 0, 0]


def inp_pri_date(start_date_str, stop_date_str):
    """(START_YEAR, START_DAY, PRISTART, PRISTOP) for the run window."""
    start_time = parse_run_datetime(start_date_str, "start_time")
    stop_time = parse_run_datetime(stop_date_str, "stop_time")

    START_YEAR = start_time.year
    START_DAY = start_time.timetuple().tm_yday

    PRISTART = [start_time.timetuple().tm_yday, start_time.hour, start_time.minute, start_time.second]

    # PRISTOP's day counts on from START_YEAR without wrapping (2023-01-01 from 2022 is day 366).
    PRISTOP_day = START_DAY + (stop_time.date() - start_time.date()).days
    PRISTOP = [PRISTOP_day, stop_time.hour, stop_time.minute, stop_time.second]

    return START_YEAR, START_DAY, PRISTART, PRISTOP

def inp_mxhist(start_time, stop_time, x_hist, mxhist_warn, segment = None):
    """(MXHIST, hint message) for histories at cadence x_hist: one day's worth, or one segment's."""
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    total_duration = stop - start
    total_seconds = total_duration.total_seconds()
    
    seconds_in_day = 86400
    seconds_in_hour = 3600
    seconds_in_min = 60
    
    n_day, n_hour, n_min, n_sec = x_hist
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec
    
    if step_seconds == 0:
        raise ValueError("Invalid prihist: step cannot be 0.")
    
    mxhist_day = seconds_in_day / step_seconds
    mxhist_hour = seconds_in_hour / step_seconds
    mxhist_min = seconds_in_min / step_seconds    
    if segment == None:
        if mxhist_day >= 1:
            mxhist_warn = (mxhist_warn + "\n" if mxhist_warn is not None else "") + f"For a Daily output set MXHIST to {int(mxhist_day)}"
        if mxhist_hour >= 1:
            mxhist_warn = mxhist_warn +  f"\nFor a Hourly output set MXHIST to {int(mxhist_hour)}"
        if mxhist_min >= 1:
            mxhist_warn = mxhist_warn + f"\nFor a Minutely output set MXHIST to {int(mxhist_min)}"
        MXHIST = mxhist_day
    else:
        segment_seconds = (segment[0] * 86400) + (segment[1] * 3600) + (segment[2] * 60) + segment[3]
        MXHIST = segment_seconds/step_seconds
        mxhist_warn = f"MXHIST minimum = 1, maximum = {int(MXHIST)} for segment run."
    return(max(1, int(MXHIST)), mxhist_warn)

def inp_sechist(SECSTART, SECSTOP, segment = None, step = 60):
    """Default SECHIST: 1 day for a window of 7 days or more, else 1 hour, fitted to the window
    and to the segment."""
    PRISTART_DAY = SECSTART[0]
    PRISTOP_DAY = SECSTOP[0]
    n_split_day = int(PRISTOP_DAY - PRISTART_DAY)
    if n_split_day >= 7:
        SECHIST = [1, 0, 0, 0]
    else:
        SECHIST = [0, 1, 0, 0]

    SECHIST = fit_cadence(SECHIST, _mtime_seconds(SECSTOP) - _mtime_seconds(SECSTART), step)
    if segment != None:
        SECHIST = fit_cadence(SECHIST, _mtime_seconds(segment), step)

    return SECHIST

def inp_prihist(PRISTART, PRISTOP, segment = None, step = 60):
    """Default PRIHIST: 1 day for 7 days or more, else 1 hour, fitted to the window; a segmented
    run uses one unit of the segment's finest field."""
    PRISTART_DAY = PRISTART[0]
    PRISTOP_DAY = PRISTOP[0]
    n_split_day = int(PRISTOP_DAY - PRISTART_DAY)
    if n_split_day >= 7:
        PRIHIST = [1, 0, 0, 0]
    else:
        PRIHIST = [0, 1, 0, 0]

    if segment != None:
        PRIHIST = _segment_cadence(segment)
    else:
        PRIHIST = fit_cadence(PRIHIST, _mtime_seconds(PRISTOP) - _mtime_seconds(PRISTART), step)
    return PRIHIST

# mkhvols.F reads each path of the OUTPUT/SECOUT 'to' form into character(len=240).
HIST_TO_PATH_MAX = 240


def _check_to_form(key, first, last):
    long_paths = [(n, p) for n, p in (("first", first), ("last", last)) if len(p) > HIST_TO_PATH_MAX]
    if long_paths:
        name, path = long_paths[0]
        raise ValueError(f"{key}: the {name} history path is {len(path)} characters, over the "
                         f"{HIST_TO_PATH_MAX} TIE-GCM can read in the '<first>','to','<last>','by','1' "
                         f"form: {path}. Use a shorter history directory (model.data.histdir) or job name.")


def inp_pri_out(start_time, stop_time, PRIHIST, MXHIST_PRIM, pri_files, histdir, run_name):
    """(OUTPUT namelist value, last primary file number) for the window; histdir is relative to
    the job's workdir or absolute."""
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')

    total_seconds = (stop - start).total_seconds()

    n_day, n_hour, n_min, n_sec = PRIHIST
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec

    data_per_file_seconds = step_seconds * int(MXHIST_PRIM)
    if data_per_file_seconds == 0:
        raise ValueError(f"Cannot size primary-history files: PRIHIST step ({step_seconds}s) * MXHIST_PRIM ({MXHIST_PRIM}) is zero.")

    number_of_files = ceil(total_seconds / data_per_file_seconds)
    if number_of_files == 1:
        pri_files_n = pri_files + 1
        OUTPUT = f"'{histdir}/{run_name}_temp_{'{:02d}'.format(pri_files)}.nc' , '{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files+1)}.nc'"
    else:
        # A continuation segment starts at pri_files + 1: prim_<pri_files> is its SOURCE and
        # output.F opens OUTPUT(1) with 'REPLACE'.
        first = pri_files + 1 if pri_files > 0 else 0
        pri_files_n = first + number_of_files
        PRIM_0 = f"{histdir}/{run_name}_prim_{'{:02d}'.format(first)}.nc"
        PRIM_N = f"{histdir}/{run_name}_prim_{'{:02d}'.format(pri_files_n)}.nc"
        _check_to_form("OUTPUT", PRIM_0, PRIM_N)
        OUTPUT = f"'{PRIM_0}','to','{PRIM_N}','by','1'"
    return OUTPUT, pri_files_n

def inp_sec_out(start_time, stop_time, SECHIST, MXHIST_SECH, sec_files, histdir, run_name):
    """(SECOUT namelist value, last secondary file number) for the window."""
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    sechist_delta = timedelta(days=SECHIST[0], hours=SECHIST[1], minutes=SECHIST[2], seconds=SECHIST[3])
    start = start + sechist_delta

    total_seconds = (stop - start).total_seconds()

    n_day, n_hour, n_min, n_sec = SECHIST
    step_seconds = (n_day * 86400) + (n_hour * 3600) + (n_min * 60) + n_sec

    data_per_file_seconds = step_seconds * int(MXHIST_SECH)
    if data_per_file_seconds == 0:
        raise ValueError(f"Cannot size secondary-history files: SECHIST step ({step_seconds}s) * MXHIST_SECH ({MXHIST_SECH}) is zero.")

    number_of_files = ceil(total_seconds / data_per_file_seconds)
    # numfiles.F counts SECSTART..SECSTOP inclusive, so exactly MXHIST_SECH + 1 histories need 2 files.
    overflow = total_seconds // step_seconds + 1 > int(MXHIST_SECH)
    sec_files_start = sec_files + 1
    sec_files_end = sec_files_start + number_of_files
    if (number_of_files == 1 or number_of_files == 0) and not overflow:
        SECOUT = f"'{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_start)}.nc'"
        sec_files_end = sec_files_start
    else:
        SECH_0 = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_start)}.nc"
        SECH_N = f"{histdir}/{run_name}_sech_{'{:02d}'.format(sec_files_end)}.nc"
        _check_to_form("SECOUT", SECH_0, SECH_N)
        SECOUT = f"'{SECH_0}','to','{SECH_N}','by','1'"

    return SECOUT, sec_files_end

def inp_sec_date(start_time, stop_time, SECHIST):
    """(SECSTART, SECSTOP): start + SECHIST to stop, days counted on from the start year."""
    start = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S')
    stop = datetime.strptime(stop_time, '%Y-%m-%dT%H:%M:%S')
    sechist_delta = timedelta(days=SECHIST[0], hours=SECHIST[1], minutes=SECHIST[2], seconds=SECHIST[3])
    sec_start = start + sechist_delta
    SECSTART_day = start.timetuple().tm_yday + (sec_start.date() - start.date()).days
    SECSTART = [SECSTART_day,sec_start.hour,sec_start.minute,sec_start.second]
    SECSTOP_day = start.timetuple().tm_yday + (stop.date() - start.date()).days
    SECSTOP = [SECSTOP_day,stop.hour,stop.minute,stop.second]

    return SECSTART, SECSTOP


def _nl_value_set(value):
    if value is None:
        return False
    if isinstance(value, (list, tuple)):
        return any(_nl_value_set(v) for v in value)
    return str(value).strip().lower() not in ("", "none", "null", "[none]", "[null]")


_IMF_PARAMS = ("BXIMF", "BYIMF", "BZIMF", "SWDEN", "SWVEL")


def _other_input_text(inp):
    return " ".join(str(x) for x in (inp.get("other_input") or []) if x is not None)


def _given(inp, name, other=None):
    """True when name is set, or given as a <name>_TIME series in other_input."""
    other = _other_input_text(inp) if other is None else other
    return _nl_value_set(inp.get(name)) or re.search(
        rf"\b{name}_TIME\b", other, re.IGNORECASE) is not None


def _potential_model(inp):
    return str(inp.get("POTENTIAL_MODEL") or "HEELIS").strip().upper()


def forcing_missing(inp):
    """[(FIELD, reason)] for the forcing input.F inp_solar requires but inp leaves unset."""
    other = _other_input_text(inp)
    pot = _potential_model(inp)
    gpi = _nl_value_set(inp.get("GPI_NCFILE"))
    missing = []
    if pot.startswith("HEELIS"):
        if gpi:
            return []
        kp = _given(inp, "KP", other)
        if not (_given(inp, "POWER", other) or kp):
            missing.append(("POWER", "POWER (or KP) is required for a HEELIS run without a GPI "
                                     "file"))
        if not (_given(inp, "CTPOTEN", other) or kp):
            missing.append(("CTPOTEN", "CTPOTEN (or KP) is required for a HEELIS run without a GPI "
                                       "file"))
        why = "for a HEELIS run without a GPI file"
        if not _given(inp, "F107", other):
            missing.append(("F107", "F107 is required " + why.format("F107")))
        if not _given(inp, "F107A", other):
            missing.append(("F107A", "F107A is required " + why.format("F107a")))
    elif pot.startswith("WEIMER"):
        if not _nl_value_set(inp.get("IMF_NCFILE")):
            for name in _IMF_PARAMS:
                if not _given(inp, name, other):
                    missing.append((name, f"{name} is required for a WEIMER run without an IMF file "
                                          f"(IMF_NCFILE): all of BXIMF, BYIMF, BZIMF, SWDEN and "
                                          f"SWVEL"))
        if not gpi:
            # input.F checks only the scalar f107/f107a here; an F107_TIME series does not count.
            why = "for a WEIMER run without a GPI file"
            if not _nl_value_set(inp.get("F107")):
                missing.append(("F107", "F107 is required " + why))
            if not _nl_value_set(inp.get("F107A")):
                missing.append(("F107A", "F107A is required " + why))
    return missing


def forcing_refusal(missing, where="the missing values"):
    """The error message for forcing_missing's list."""
    names = [f for f, _ in missing]
    files = []
    if any(f in names for f in ("POWER", "CTPOTEN", "F107", "F107A")):
        files.append("a GPI_NCFILE (or 'gen')")
    if any(f in names for f in _IMF_PARAMS):
        files.append("an IMF_NCFILE (or 'gen')")
    return ("Incomplete forcing, not writing a namelist TIE-GCM would reject: "
            + "; ".join(why for _, why in missing)
            + f". Provide {' and '.join(files) or 'the forcing files'}, or {where}.")


def forcing_conflicts(inp):
    """[(FIELD, reason)] for forcing values inp gives that input.F inp_solar refuses."""
    other = _other_input_text(inp)
    if _potential_model(inp).startswith("HEELIS"):
        if _nl_value_set(inp.get("IMF_NCFILE")):
            return [("IMF_NCFILE", "IMF_NCFILE is read only by a WEIMER run: set inp.IMF_NCFILE to "
                                   "null for HEELIS")]
        return []
    if not _potential_model(inp).startswith("WEIMER"):
        return []
    conflicts = []
    kp = _given(inp, "KP", other)
    if kp:
        conflicts.append(("KP", "KP cannot be given for a WEIMER run (the model stops): set inp.KP "
                                "to null"))
    # input.F: `ctpoten /= spval .or. ntimes_ctpoten > 0 .and. kp == spval ...` (.and. first)
    if _nl_value_set(inp.get("CTPOTEN")) or (_given(inp, "CTPOTEN", other) and not kp):
        conflicts.append(("CTPOTEN", "CTPOTEN cannot be given for a WEIMER run (the Weimer model "
                                     "computes it): set inp.CTPOTEN to null"))
    if _nl_value_set(inp.get("IMF_NCFILE")) and all(_given(inp, n, other) for n in _IMF_PARAMS):
        conflicts.append(("IMF_NCFILE", "an IMF run cannot also give all of BXIMF, BYIMF, BZIMF, "
                                        "SWDEN and SWVEL (at least one must come from the IMF "
                                        "file): set one of them to null"))
    return conflicts


# input.F requires SOURCE_START's time of day to equal PRISTART's, and nc_rdhist needs
# SOURCE_START to match a history of SOURCE exactly.
def hms_str(mtime):
    return ":".join(f"{int(x):02d}" for x in list(mtime)[1:4])


def mtime_ints(value):
    """A 'd h m s' / [d, h, m, s] model time as a list of ints, or None."""
    try:
        if isinstance(value, (list, tuple)):
            return [int(x) for x in value]
        return [int(x) for x in str(value).replace("'", "").replace(",", " ").split()]
    except (TypeError, ValueError):
        return None


def hms_differs(source_start, pristart):
    ss = mtime_ints(source_start)
    return ss is not None and len(ss) == 4 and ss[1:4] != [int(x) for x in list(pristart)[1:4]]


def source_start_hms_error(source_start, pristart):
    if not hms_differs(source_start, pristart):
        return None
    return (f"SOURCE_START {source_start} must have the run start's time of day "
            f"{hms_str(pristart)} (TIE-GCM requires it; only the day may differ).")


def source_time_of_day_error(source, source_mtimes, pristart, start_time=None):
    """(prompt message, replay refusal) when SOURCE has no history at PRISTART's time of day,
    else None."""
    if any(list(m[1:4]) == list(pristart[1:4]) for m in source_mtimes):
        return None
    hms = hms_str(pristart)
    have = sorted({hms_str(m) for m in source_mtimes})
    day = f"{start_time[:11]}" if start_time else ""
    advice = (f"Re-run with a start time at one of those times of day (e.g. "
              f"{day}{have[0] if have else '00:00:00'}), or give a SOURCE "
              f"file with a history at {hms}.")
    message = (f"SOURCE {source} has no history at the run start's time of day {hms} "
               f"(its histories are at {', '.join(have) if have else 'no time'}). TIE-GCM needs a "
               f"start history at the run start's time of day (only the day may differ). "
               f"{advice}")
    refusal = (f"No SOURCE history at the run start's time of day {hms} (SOURCE {source} "
               f"has {', '.join(have) if have else 'none'}); TIE-GCM would stop at start-up. "
               f"{advice}")
    return message, refusal


def resolve_source_start(source, source_start, source_mtimes, what="inp.SOURCE_START", pristart=None):
    """(SOURCE_START, note or None), checked against the histories of SOURCE.

    Unset or not a history: derived from SOURCE's only (or first at PRISTART's time of day)
    history; ValueError when it is not a history and SOURCE holds several."""
    histories = [[int(x) for x in m] for m in source_mtimes]
    listed = ", ".join(" ".join(map(str, m)) for m in histories) or "none"
    if not histories:
        raise ValueError(f"SOURCE {source} holds no history (no mtime): TIE-GCM cannot start from it.")
    if source_start is None or str(source_start).strip().lower() in ("", "none", "null"):
        at_start = [m for m in histories if pristart is not None and m[1:4] == [int(x) for x in list(pristart)[1:4]]]
        return " ".join(map(str, (at_start or histories)[0])), None
    ss = mtime_ints(source_start)
    if ss is not None and ss in histories:
        return source_start, None
    if len(histories) == 1:
        value = " ".join(map(str, histories[0]))
        return value, (f"NOTE: {what} {source_start} is not a history of SOURCE {source}; SOURCE "
                       f"holds one history, so {what} is re-derived as {value}.")
    raise ValueError(
        f"{what} {source_start} is not a history of SOURCE {source} (TIE-GCM would stop with "
        f"'Source history not found'): SOURCE holds several histories ({listed}); set "
        f"{what} to the one to start from.")


def date_order_problem(label, value, after=None, after_label=None, strictly_after=True,
                       before=None, before_label=None, strictly_before=False):
    """An error message when value is not after `after` / before `before`, else None."""
    try:
        t = parse_run_datetime(value, label)
        if after is not None:
            a = parse_run_datetime(after, after_label)
            if (t <= a) if strictly_after else (t < a):
                return (f"{label} {value} must be {'after' if strictly_after else 'on or after'} "
                        f"{after_label} {after}.")
        if before is not None:
            b = parse_run_datetime(before, before_label)
            if (t >= b) if strictly_before else (t > b):
                return (f"{label} {value} must be {'before' if strictly_before else 'on or before'} "
                        f"{before_label} {before}.")
    except ValueError as e:
        return str(e)
    return None


def mxday(year):
    """The last model day a run starting in `year` may reach (input.F: 367 leap, else 366)."""
    return 367 if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0 else 366


def year_boundary(start_time, stop_time):
    """(start year, PRISTOP day, mxday) when the run passes the start year's last model day."""
    year, _, _, pristop = inp_pri_date(start_time, stop_time)
    if pristop[0] > mxday(year):
        return year, pristop[0], mxday(year)
    return None


def window_problems(inp):
    """Problems with the run and secondary windows and the segment."""
    problems = []
    seg_value = inp.get("segment")
    segment = None
    if _nl_value_set(seg_value):
        try:
            segment = parse_segment(seg_value)
            if segment is None:
                problems.append(f"inp.segment {seg_value!r} is all zero: give a segment length or none")
        except ValueError as e:
            problems.append(f"inp.{e}")
    start, stop = inp.get("start_time"), inp.get("stop_time")
    if not (_nl_value_set(start) and _nl_value_set(stop)):
        return problems                       # benchmark: the window is in the namelist keys
    err = date_order_problem("inp.stop_time", stop, start, "inp.start_time")
    if err:
        return problems + [err]
    sec_start, sec_stop = inp.get("secondary_start_time"), inp.get("secondary_stop_time")
    if _nl_value_set(sec_start):
        err = date_order_problem("inp.secondary_start_time", sec_start, start, "inp.start_time",
                                 strictly_after=False, before=stop, before_label="inp.stop_time",
                                 strictly_before=True)
        if err:
            problems.append(err)
        elif _nl_value_set(sec_stop):
            err = date_order_problem("inp.secondary_stop_time", sec_stop, sec_start,
                                     "inp.secondary_start_time", before=stop,
                                     before_label="inp.stop_time")
            if err:
                problems.append(err)
    crossing = year_boundary(start, stop)
    if crossing and segment is None and not problems:
        year, day, last = crossing
        problems.append(f"inp.segment is required: the run crosses the {year}->{year + 1} year "
                        f"boundary (stop day {day} > {last})")
    return problems


LBC_OTHER_KEYS = ("CTMT_NCFILE", "SABER_NCFILE", "TIDI_NCFILE", "TIDE", "TIDE2")
_GSWM_KEYS = ("GSWM_MI_DI_NCFILE", "GSWM_MI_SDI_NCFILE", "GSWM_NM_DI_NCFILE", "GSWM_NM_SDI_NCFILE")
_NUDGE_KEYS = ("NUDGE_NCPRE", "NUDGE_NCPOST", "NUDGE_FLDS", "NUDGE_LBC", "NUDGE_F4D",
               "NUDGE_USE_REFDATE", "NUDGE_REFDATE", "NUDGE_SPONGE", "NUDGE_DELTA", "NUDGE_POWER",
               "NUDGE_ALPHA")


def _tide_set(value):
    """True for TIDE / TIDE2 with a non-zero amplitude."""
    if not _nl_value_set(value):
        return False
    items = value if isinstance(value, (list, tuple)) else [value]
    for item in items:
        for x in re.split(r"[,\s]+", str(item).strip().strip("[]")):
            x = x.strip("'\"")
            if not x:
                continue
            try:
                if float(x.lower().replace("d", "e")) != 0:
                    return True
            except ValueError:
                return True
    return False


def _lbc_families(inp):
    """{family: [keys set]} of the lower-boundary tide sources inp gives."""
    fams = {"GSWM": [k for k in _GSWM_KEYS if _nl_value_set(inp.get(k))],
            "TIDE": [k for k in ("TIDE", "TIDE2") if _tide_set(inp.get(k))],
            "SABER/TIDI": [k for k in ("SABER_NCFILE", "TIDI_NCFILE") if _nl_value_set(inp.get(k))],
            "CTMT": [k for k in ("CTMT_NCFILE",) if _nl_value_set(inp.get(k))]}
    return {f: keys for f, keys in fams.items() if keys}


def lbc_other_set(inp):
    return any(f != "GSWM" for f in _lbc_families(inp))


def lbc_conflicts(inp, horires):
    """[(KEY, reason)] for lower-boundary and nudging input the model refuses."""
    out = []
    fams = _lbc_families(inp)
    if len(fams) > 1:
        keys = [k for f in fams.values() for k in f]
        out.append((keys[0], f"inp.{', inp.'.join(keys)} give {len(fams)} lower-boundary tide "
                             f"sources ({', '.join(fams)}): keep one and set the others to null"))
    try:
        five = float(horires) == 5.0
    except (TypeError, ValueError):
        five = False
    for key in ("BGRDDATA_NCFILE", "CTMT_NCFILE"):
        if five and _nl_value_set(inp.get(key)):
            out.append((key, f"inp.{key} is not supported on the 5-degree grid: set it to null"))
    nudge = [k for k in _NUDGE_KEYS if _nl_value_set(inp.get(k))]
    if nudge and not _nl_value_set(inp.get("NUDGE_NCFILE")):
        out.append((nudge[0], f"inp.{', inp.'.join(nudge)} set without inp.NUDGE_NCFILE"))
    return out


def _cadence_value(value):
    ints = mtime_ints(value)
    return ints if ints is not None and len(ints) == 4 and any(ints) else None


def cadence_problems(inp, segment, step):
    """Problems with STEP and with PRIHIST / SECHIST against STEP and their window or segment."""
    problems = []
    try:
        step = int(float(step))
    except (TypeError, ValueError):
        return [f"inp.STEP {step!r} is not a whole number of seconds"]
    if step <= 0:
        return [f"inp.STEP {step} must be > 0"]
    seg_s = _mtime_seconds(segment) if segment else None
    windows = (("PRIHIST", "PRISTART", "PRISTOP"), ("SECHIST", "SECSTART", "SECSTOP"))
    for key, w0, w1 in windows:
        cadence = _cadence_value(inp.get(key))
        if cadence is None:
            continue
        if seg_s:
            span, where = seg_s, f"the {' '.join(map(str, segment))} segment"
        else:
            a, b = _cadence_value(inp.get(w0)), _cadence_value(inp.get(w1))
            if a is None or b is None:
                continue
            span, where = _mtime_seconds(b) - _mtime_seconds(a), f"{w0}..{w1}"
        try:
            fitted = fit_cadence(cadence, span, step)
        except ValueError as e:
            problems.append(f"inp.{key}: {e}")
            continue
        if list(fitted) != list(cadence):
            problems.append(f"inp.{key} {' '.join(map(str, cadence))} does not fit STEP {step} s and "
                            f"{where}: use {' '.join(map(str, fitted))}")
    return problems
