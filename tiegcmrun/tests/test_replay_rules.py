"""-o / engage -to replay rules: the checks the prompts make hold for a replayed options JSON too."""
import argparse
import json
import os
import socket
import sys
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

import tiegcmrun  # noqa: E402
import misc  # noqa: E402
import engage_solver  # noqa: E402
import replay_solver  # noqa: E402
from namelist_solver import forcing_conflicts  # noqa: E402
from output_solver import apply_system_config  # noqa: E402

START, STOP = "2016-08-02T00:00:00", "2016-08-09T11:00:00"
needs_env = pytest.mark.skipif(not (os.environ.get("TIEGCMHOME") and os.environ.get("TIEGCMDATA")),
                               reason="needs TIEGCMHOME and TIEGCMDATA")


@pytest.fixture
def od():
    with open(tiegcmrun.OPTION_DESCRIPTIONS_FILE, encoding="utf-8") as f:
        return apply_system_config(json.load(f), socket.gethostname())


def options_json(tmp_path, **inp):
    """A small standalone aitken options JSON in tmp_path; the run directory is tmp_path/run."""
    options = {"simulation": {"job_name": "x", "hpc_system": "aitken"},
               "model": {"data": {"parentdir": "run"}, "specification": {"horires": 2.5}},
               "inp": {"start_time": START, "stop_time": STOP, **inp},
               "job": {"project_code": "s1234", "queue": "normal"}}
    path = tmp_path / "t.json"
    path.write_text(json.dumps(options))
    return path


def run_main(path, *flags):
    return tiegcmrun.main(["-o", str(path), "--no-tui", *flags])


# --- B02: dated GPI / IMF files must cover the run --------------------------------------------------

def test_validate_options_refuses_stale_dated_gpi_and_imf(tmp_path):
    inp = {"start_time": START, "stop_time": STOP, "GPI_NCFILE": "/d/gpi_2010001-2010365.nc",
           "IMF_NCFILE": "/d/imf_OMNI_2017001-2017365.nc"}
    problems = replay_solver._cover_problems(inp)
    assert problems == [
        "inp.GPI_NCFILE gpi_2010001-2010365.nc covers 2010001..2010365, not this run "
        "(2016-08-02T00:00:00..2016-08-09T11:00:00): set it to 'gen' (a file is generated for these "
        "dates) or to a file that covers them.",
        "inp.IMF_NCFILE imf_OMNI_2017001-2017365.nc covers 2017001..2017365, not this run "
        "(2016-08-02T00:00:00..2016-08-09T11:00:00): set it to 'gen' (a file is generated for these "
        "dates) or to a file that covers them."]
    # 'gen', a $-value, an undated name, unset, a generated key and a benchmark JSON pass
    for value in ("gen", "$TIEGCMDATA/gpi_2010001-2010365.nc", "/d/gpi_test.nc", None):
        assert replay_solver._cover_problems(dict(inp, GPI_NCFILE=value, IMF_NCFILE=None)) == []
    assert replay_solver._cover_problems(inp, ("inp.GPI_NCFILE", "inp.IMF_NCFILE")) == []
    assert replay_solver._cover_problems({"GPI_NCFILE": "/d/gpi_2010001-2010365.nc"}) == []
    # the prompts keep their wording
    assert misc.forcing_cover_problem("GPI_NCFILE", "gpi_2010001-2010365.nc", START, STOP).endswith(
        "enter 'gen' to generate a file for these dates, or a file that covers them.")


@needs_env
def test_replay_refuses_stale_gpi_end_to_end(tmp_path, capsys):
    stale = tmp_path / "gpi_2010001-2010365.nc"
    stale.write_text("")
    assert run_main(options_json(tmp_path, GPI_NCFILE=str(stale)), "--check") == 1
    out = capsys.readouterr().out
    assert "inp.GPI_NCFILE gpi_2010001-2010365.nc covers 2010001..2010365, not this run" in out
    assert not (tmp_path / "run").exists()


# --- B03: KP with both POWER and CTPOTEN ------------------------------------------------------------

def test_forcing_conflicts_kp_power_ctpoten():
    manual = {"KP": 3, "POWER": 20, "CTPOTEN": 40, "F107": 70, "F107A": 70}
    for inp in (manual, dict(manual, GPI_NCFILE="/d/gpi.nc"), dict(manual, POTENTIAL_MODEL="HEELIS")):
        assert forcing_conflicts(inp) == [
            ("KP", "KP cannot be given together with both POWER and CTPOTEN (the model stops: KP stands in "
                   "for the one that is missing): set inp.KP, or one of inp.POWER and inp.CTPOTEN, to null")]
    assert forcing_conflicts({"KP": 3, "POWER": 20}) == []
    assert forcing_conflicts({"GPI_NCFILE": "/d/gpi.nc", "KP": 3}) == []
    # a CTPOTEN_TIME series counts as given
    assert forcing_conflicts({"KP": 3, "POWER": 20, "other_input": ["CTPOTEN_TIME = 1,2"]})[0][0] == "KP"
    # WEIMER keeps its own rules
    assert [k for k, _ in forcing_conflicts(dict(manual, POTENTIAL_MODEL="WEIMER"))] == ["KP", "CTPOTEN"]


def test_imf_conflict_text_starts_with_key():
    inp = {"POTENTIAL_MODEL": "WEIMER", "IMF_NCFILE": "/d/imf.nc", "BXIMF": 0, "BYIMF": 0, "BZIMF": -5,
           "SWDEN": 5, "SWVEL": 400, "GPI_NCFILE": "/d/gpi.nc"}
    problems = replay_solver._forcing_problems(inp)
    assert problems == ["inp.IMF_NCFILE is set: an IMF run cannot also give all of BXIMF, BYIMF, BZIMF, "
                        "SWDEN and SWVEL (at least one must come from the IMF file): set one of them to null"]


def test_prompt_loop_drops_prefilled_forcing_under_gpi(capsys):
    unused = []
    prefilled = {"default": 3, "_prefilled": True}
    assert tiegcmrun.skipped_forcing_value("KP", prefilled, None, unused) is None
    assert tiegcmrun.skipped_forcing_value("POWER", {"default": None}, None, unused) is None
    assert unused == ["inp.KP 3"]
    # a benchmark keeps its values
    assert tiegcmrun.skipped_forcing_value("KP", prefilled, "nov2003_heelis_gpi", unused) == 3
    assert unused == ["inp.KP 3"]
    assert tiegcmrun._forcing_source({"GPI_NCFILE": "/d/gpi_x.nc", "KP": 3}) == "GPI_NCFILE gpi_x.nc"
    assert tiegcmrun._forcing_source({"KP": 3}) == "KP 3"
    assert tiegcmrun._forcing_source({"POTENTIAL_MODEL": "WEIMER"}) == "POTENTIAL_MODEL WEIMER"
    assert (tiegcmrun.unused_forcing_note(["inp.KP 3", "inp.POWER 20"], {"GPI_NCFILE": "/d/gpi_x.nc"})
            == "NOTE: inp.KP 3, inp.POWER 20 (prefilled) not used: the forcing comes from GPI_NCFILE gpi_x.nc.")


# --- B04 / B30: --check refuses a pending 'gen' without gcmprocpy; rules before generation ----------

@needs_env
def test_check_refuses_gen_without_gcmprocpy(tmp_path, monkeypatch, capsys):
    path = options_json(tmp_path, GPI_NCFILE="gen")
    monkeypatch.setitem(sys.modules, "gcmprocpy", None)          # import fails
    assert run_main(path, "--check") == 1
    out = capsys.readouterr().out
    assert "ERROR: inp.GPI_NCFILE = 'gen': GPI generation needs the 'gcmprocpy' package" in out
    assert "Nothing was written: install gcmprocpy, or set inp.GPI_NCFILE to a file that covers the run." in out
    assert not (tmp_path / "run").exists()
    monkeypatch.setitem(sys.modules, "gcmprocpy", types.ModuleType("gcmprocpy"))
    assert run_main(path, "--check") == 0
    assert "--check: inp.GPI_NCFILE is 'gen': not generated now" in capsys.readouterr().out
    assert not (tmp_path / "run").exists()


@needs_env
def test_conflict_checked_before_generation(tmp_path, monkeypatch, capsys):
    path = options_json(tmp_path, GPI_NCFILE="gen", KP=3, POWER=20, CTPOTEN=40)
    called = []
    monkeypatch.setattr(misc, "generate_pending", lambda *a: called.append(a) or pytest.fail("generated"))
    monkeypatch.setitem(sys.modules, "gcmprocpy", types.ModuleType("gcmprocpy"))
    assert run_main(path) == 1
    assert "ERROR: inp.KP cannot be given together with both POWER and CTPOTEN" in capsys.readouterr().out
    assert called == [] and not (tmp_path / "run").exists()


def test_generate_pending_message_single_period(monkeypatch):
    def boom(*args):
        raise RuntimeError("no OMNI data for these dates.")
    monkeypatch.setattr(misc, "generate_gpi_file", boom)
    with pytest.raises(RuntimeError) as e:
        misc.generate_pending({"GPI_NCFILE": "gen", "start_time": START, "stop_time": STOP}, "/tmp/w")
    assert "generation failed: no OMNI data for these dates. Provide a file path instead" in str(e.value)
    assert ".." not in str(e.value).split("Provide")[0]


def test_gcmprocpy_problem(monkeypatch):
    monkeypatch.setitem(sys.modules, "gcmprocpy", None)
    text = misc.gcmprocpy_problem("IMF")
    assert text.startswith("IMF generation needs the 'gcmprocpy' package (>=1.5) in the active environment "
                           "(import failed: ")
    with pytest.raises(RuntimeError, match="generation needs the 'gcmprocpy' package"):
        misc._import_gcmprocpy("IMF")
    monkeypatch.setitem(sys.modules, "gcmprocpy", types.ModuleType("gcmprocpy"))
    assert misc.gcmprocpy_problem("IMF") is None


# --- B08: null sections under engage ------------------------------------------------------------------

def test_section_problem_null_under_engage():
    for key in ("inp", "job"):
        assert replay_solver.section_problem({key: None}) is None                  # standalone -o
        assert replay_solver.section_problem({key: None}, allow_null=False) == (
            f"{key} is null: it must be a JSON object (delete the key, or give {{}}, for the defaults)")
        assert tiegcmrun.options_file_problem({key: None}, allow_null=False).startswith(f"{key} is null")
    assert replay_solver.section_problem({"job": {"resource": None}}, allow_null=False).startswith(
        "job.resource is null")
    assert replay_solver.section_problem({"inp": 3}, allow_null=False) == "inp must be a JSON object"
    assert replay_solver.gpi_absent(None)


# --- B09: a missing or unreadable SOURCE under engage --------------------------------------------------

def test_engage_source_missing_one_error(tmp_path, od):
    def error(source):
        with pytest.raises(ValueError) as e:
            engage_solver.source_mtimes({"SOURCE": source}, od)
        return misc_plain(str(e.value))
    assert error("/no/such/source.nc") == "inp.SOURCE = '/no/such/source.nc' not found"
    assert error(str(tmp_path)) == f"inp.SOURCE = {str(tmp_path)!r} is a directory, not a file"
    bad = tmp_path / "bad.nc"
    bad.write_text("not netCDF")
    assert error(str(bad)).startswith(f"inp.SOURCE = {str(bad)!r} cannot be read as a TIE-GCM history file (")


def misc_plain(text):
    import re
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


@needs_env
def test_engage_source_bare_name_under_tiegcmdata(od):
    data = os.environ["TIEGCMDATA"]
    found = misc.find_file("decsol_f70.nc", data)
    if not found:
        pytest.skip("no decsol_f70.nc under TIEGCMDATA")
    inp = {"SOURCE": "decsol_f70.nc"}
    assert engage_solver.source_mtimes(inp, od)
    assert inp["SOURCE"] == found


# --- B10: bare data file names ----------------------------------------------------------------------

def test_resolve_json_paths_bare_names(tmp_path):
    (tmp_path / "beside.nc").write_text("")
    options = {"inp": {"SOURCE": "decsol_f70.nc", "GPI_NCFILE": "beside.nc", "IMF_NCFILE": "sub/imf.nc",
                       "AMIENH": "/abs/amie.nc", "HE_COEFS_NCFILE": "gen"},
               "model": {"data": {"workdir": "stdout", "modelexe": "exe"}}}
    misc.resolve_json_paths(options, str(tmp_path))
    inp = options["inp"]
    assert inp["SOURCE"] == "decsol_f70.nc"                                  # looked up by coerce
    assert inp["GPI_NCFILE"] == str(tmp_path / "beside.nc")                  # beside the JSON
    assert inp["IMF_NCFILE"] == str(tmp_path / "sub" / "imf.nc")             # a relative path
    assert inp["AMIENH"] == "/abs/amie.nc" and inp["HE_COEFS_NCFILE"] == "gen"
    assert options["model"]["data"]["workdir"] == str(tmp_path / "stdout")  # run paths as before


# --- B24 (a), B27: SECFLDS --------------------------------------------------------------------------

def test_secflds_replay_empty_list_refused():
    with pytest.raises(ValueError) as e:
        replay_solver.secflds_replay({"inp": {"SECFLDS": []}}, where="inp.SECFLDS (TIE-GCM options JSON)")
    assert str(e.value) == ("inp.SECFLDS (TIE-GCM options JSON) = []: an empty list writes no SECFLDS; give "
                            "at least one secondary-history field (the valid names are in "
                            "tiegcmrun/secflds_fields.json), or delete the key for the default list")
    for value in (None, [None]):
        options = {"inp": {"SECFLDS": value}}
        replay_solver.secflds_replay(options)
        assert options["inp"]["SECFLDS"] == value


def test_secflds_string_note_names_split(capsys):
    options = {"inp": {"SECFLDS": "TN UN VN"}}
    replay_solver.secflds_replay(options, where="inp.SECFLDS (TIE-GCM options JSON)")
    assert options["inp"]["SECFLDS"] == ["TN", "UN", "VN"]
    assert ("NOTE: inp.SECFLDS (TIE-GCM options JSON) normalised: 'TN UN VN' split into TN, UN, VN."
            in capsys.readouterr().out)
    options = {"inp": {"SECFLDS": "TN"}}
    replay_solver.secflds_replay(options, where="inp.SECFLDS")
    assert options["inp"]["SECFLDS"] == ["TN"]
    assert "NOTE: inp.SECFLDS normalised: 'TN' read as the list ['TN']." in capsys.readouterr().out


# --- B15: duplicate keys ----------------------------------------------------------------------------

def test_read_options_json_duplicate_key(tmp_path):
    path = tmp_path / "tiegcmrun_input.json"
    path.write_text('{"inp": {"KP": 3, "POWER": 20, "KP": 4}}')
    with pytest.raises(ValueError) as e:
        replay_solver.read_options_json(str(path), "tiegcmrun_input.json")
    assert str(e.value) == ("tiegcmrun_input.json: duplicate key 'KP' (JSON would keep only the last value): "
                            "keep one of them")
    path.write_text('{"inp": {"KP": 3}, "job": {"KP": 4}}')                  # other objects: fine
    assert replay_solver.read_options_json(str(path))["job"] == {"KP": 4}
    path.write_text('{"inp": ')
    with pytest.raises(ValueError, match="is not valid JSON"):
        replay_solver.read_options_json(str(path))


# --- B26: one saved-options writer ------------------------------------------------------------------

def test_replay_options_for_json():
    options = {"job": {"queue": "normal"}, "coupled_unset": ["OPFLOOR"],
               "inp": {"PRIHIST": "0 1 0 0", "MXHIST_PRIM": 1, "OUTPUT": "o", "SECHIST": "0 1 0 0",
                       "MXHIST_SECH": 1, "PRISTART": "1 0 0 0", "PRISTOP": "2 0 0 0", "SECSTART": "x",
                       "SECSTOP": "y", "SECOUT": "s", "SOURCE": "/d/src.nc", "START_YEAR": 2016},
               "simulation": {"job_name": "x"},
               "model": {"data": {"make": "m", "utildir": "u", "log_file": "l", "input_file": None,
                                  "modeldir": "/t"},
                         "specification": {"horires": 2.5, "segmentation": True}}}
    out = engage_solver.replay_options_for_json(options)
    assert list(out) == ["simulation", "model", "inp", "job", "coupled_unset"]
    assert out["inp"] == {"SOURCE": "/d/src.nc", "START_YEAR": 2016}
    assert out["model"] == {"data": {"modeldir": "/t"}, "specification": {"horires": 2.5}}
    assert options["inp"]["PRIHIST"] == "0 1 0 0"                            # a copy


# --- B24 (b), (d), (e): messages --------------------------------------------------------------------

def test_exe_problems_names_zitop():
    identity = {"path": "/b/tiegcm.exe", "exists": True, "elf": True, "nlevp1": 57, "vertres": 0.25,
                "nlon": 144, "horires": 2.5, "nmlat": None, "mage_oneway_syms": ["x"]}
    _, problems = misc.exe_problems(identity, "model.data.modelexe", True, 2.5, 0.25, 5, zitop=9)
    assert len(problems) == 1
    assert ("was built for vertres 0.25 (nlevp1 57), the run needs 0.25 (nlevp1 65) for zitop 9; at vertres "
            "0.25 this build's nlevp1 57 is zitop 7: set model.specification.zitop to that, or rebuild"
            in problems[0])
    # the default zitop: the vertres wording alone
    _, problems = misc.exe_problems(dict(identity, nlevp1=113, vertres=0.125), "model.data.modelexe", True,
                                    2.5, 0.25, 5, zitop=7)
    assert "was built for vertres 0.125 (nlevp1 113), the run needs 0.25 (nlevp1 57) - " in problems[0]


def test_exe_dollar_refused(capsys):
    options = {"model": {"data": {"modelexe": "$TIEGCMHOME/tiegcm_2.5/tiegcm.exe"},
                         "specification": {"horires": 2.5, "vertres": 0.25, "nres_grid": 5}}}
    expanded = os.path.expandvars("$TIEGCMHOME/tiegcm_2.5/tiegcm.exe")
    args = argparse.Namespace(onlycompile=False, compile=False, engage={})
    with pytest.raises(ValueError) as e:
        tiegcmrun.check_executables(options, args, "/run")
    assert ("model.data.modelexe: $TIEGCMHOME/tiegcm_2.5/tiegcm.exe not found: environment variables are "
            f"not expanded in the options JSON; give the path ({expanded})") in str(e.value)
    assert "cannot run this coupled run" in str(e.value)
    # standalone: a warning, as for a missing build
    args.engage = None
    warned = tiegcmrun.check_executables(options, args, "/run")
    assert warned and "environment variables are not expanded" in warned[0]


def test_main_catches_oserror(monkeypatch, capsys):
    def missing(argv):
        raise FileNotFoundError(2, "No such file or directory", "/d/nov2003.nc")
    monkeypatch.setattr(tiegcmrun, "tiegcmrun", missing)
    assert tiegcmrun.main(["-bench", "nov2003_heelis_gpi"]) == 1
    assert "ERROR: [Errno 2] No such file or directory: '/d/nov2003.nc'" in capsys.readouterr().out


# --- B06 (-o and engage -to): another machine's queue and job fields ------------------------------------

def test_replay_refuses_foreign_queue(od, capsys):
    options = {"simulation": {"job_name": "x", "hpc_system": "aitken"},
               "job": {"queue": "main", "project_code": "s1234"}}
    with pytest.raises(ValueError) as e:
        replay_solver.settle_machine_job(options, od)
    text = misc_plain(str(e.value))
    assert text.startswith("job.queue='main' is not a queue of aitken (config/machines.yaml queues: ")
    assert text.endswith("set job.queue to one of them, or pass --rederive job to take the machine's default.")
    options["job"] = {"queue": "normal", "job_priority": "economy", "project_code": "s1234"}
    notes = replay_solver.settle_machine_job(options, od)
    assert notes == ["NOTE: job.job_priority='economy' dropped: not a job field of aitken (config/machines.yaml)."]
    assert options["job"] == {"queue": "normal", "project_code": "s1234"}
    assert notes[0] in capsys.readouterr().out
    # --rederive job drops the block first: nothing to refuse
    options = {"simulation": {"hpc_system": "aitken"}, "job": {"queue": "main"}}
    replay_solver.rederive_pop(options, ["job"])
    assert replay_solver.settle_machine_job(options, od) == []


# --- B22: aitken's default node type ----------------------------------------------------------------

def test_machines_yaml_node_type_default_rome():
    assert misc.machines_config()["aitken"]["node_type_default"] == "rome"


def test_engage_drops_foreign_job_field(od, monkeypatch, tmp_path, capsys):
    """engage -to: a job field of another machine is dropped with the -o NOTE, not hinted."""
    monkeypatch.setattr(engage_solver, "source_mtimes", lambda inp, od: [[214, 0, 0, 0]])
    monkeypatch.setattr(engage_solver, "resolve_source_start", lambda *a, **k: ([214, 0, 0, 0], None))
    monkeypatch.setattr(engage_solver, "replay_exes", lambda data, engage: None)
    options = {"simulation": {"job_name": "x", "hpc_system": "aitken"},
               "model": {"data": {}, "specification": {}},
               "inp": {"SOURCE": "/d/source.nc"},
               "job": {"queue": "normal", "job_priority": "economy"}}
    engage = {"job_name": "x", "hpc_system": "aitken", "parentdir": str(tmp_path), "horires": 2.5,
              "horires_coupled": 2.5, "skip": [], "start_time": START, "stop_time": STOP,
              "segment": [1, 0, 0, 0], "queue": "normal", "walltime": "01:00:00", "model": "rom_ait",
              "project_code": "s1234"}
    engage_solver.engage_options_updater(options, engage, od)
    out = misc_plain(capsys.readouterr().out)
    assert "job_priority" not in options["job"]
    assert "NOTE: job.job_priority='economy' dropped: not a job field of aitken (config/machines.yaml)." in out
    assert "but job has" not in out


# --- the verifier's regression: a file this run generated is not cover-checked ----------------------

def test_generated_forcing_file_passes_cover_check(tmp_path, monkeypatch):
    lagging = str(tmp_path / "gpi_2016214-2016218.nc")            # data stop short of STOP (2016222)
    inp = {"start_time": START, "stop_time": STOP, "GPI_NCFILE": lagging}
    monkeypatch.setattr(misc, "GENERATED_FORCING_FILES", set())
    assert len(replay_solver._cover_problems(inp)) == 1            # not generated here: refused
    misc.GENERATED_FORCING_FILES.add(os.path.abspath(lagging))
    assert replay_solver._cover_problems(inp) == []
    # end to end through generate_pending: the generation site warns, validation does not refuse
    monkeypatch.setattr(misc, "GENERATED_FORCING_FILES", set())

    class FakeGcmprocpy:
        @staticmethod
        def generate_gpi(**kw):
            return None

        @staticmethod
        def save_gpi(ds, output_dir):
            path = os.path.join(output_dir, "gpi_2016214-2016218.nc")
            open(path, "w").close()
            return path
    monkeypatch.setattr(misc, "_import_gcmprocpy", lambda kind: FakeGcmprocpy)
    inp = {"start_time": START, "stop_time": STOP, "GPI_NCFILE": "gen"}
    misc.generate_pending(inp, str(tmp_path))
    assert inp["GPI_NCFILE"] == lagging
    assert replay_solver._cover_problems(inp) == []
