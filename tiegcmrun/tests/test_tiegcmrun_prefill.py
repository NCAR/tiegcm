"""--prefill: the prompts / wizard start from a saved tiegcmrun_parameters.json's answers."""
import argparse
import builtins
import copy
import json
import os
import socket
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

import tiegcmrun  # noqa: E402
from output_solver import apply_system_config  # noqa: E402


def saved_parameters(data_dir):
    """A tiegcmrun_parameters.json as a standalone derecho run saves it in <run>/stdout."""
    gpi = os.path.join(data_dir, "gpi_test.nc")
    open(gpi, "w").close()
    return {
        "simulation": {"job_name": "storm", "hpc_system": "derecho"},
        "model": {
            "data": {"modeldir": "/old/tiegcm", "parentdir": "..", "execdir": "../exec", "workdir": ".",
                     "histdir": "../hist", "utildir": "/old/tiegcm/scripts", "tgcmdata": "/old/data",
                     "input_file": "./storm_1.25x0.125.inp", "log_file": "./storm.out",
                     "make": "/old/tiegcm/scripts/Make.intel_de", "modelexe": "../exec/tiegcm.exe"},
            "specification": {"horires": 1.25, "vertres": "0.125", "zitop": 7, "mres": 1,
                              "nres_grid": 6, "segmentation": False}},
        "inp": {"LABEL": "storm_1.25x0.125", "start_time": "2015-03-16T00:00:00",
                "stop_time": "2015-03-19T00:00:00", "secondary_start_time": "2015-03-16T00:00:00",
                "secondary_stop_time": "2015-03-19T00:00:00", "solar_flux_level": "high",
                "SOURCE": "/old/data/source_start/mareqx_f200.nc", "segment": [None],
                "SOURCE_START": "81 0 0 0", "START_YEAR": 2015, "START_DAY": 75,
                "CALENDAR_ADVANCE": 1, "PRISTART": "75 0 0 0", "PRISTOP": "78 0 0 0",
                "PRIHIST": "0 1 0 0", "MXHIST_PRIM": 24, "OUTPUT": "'../hist/x_prim_01.nc'",
                "STEP": 10, "NSTEP_SUB": 10, "SECHIST": "0 0 30 0", "MXHIST_SECH": 48,
                "SECSTART": "75 0 30 0", "SECSTOP": "78 0 0 0", "SECOUT": "'../hist/x_sech_01.nc'",
                "SECFLDS": ["TN", "UN", "VN"], "ELECTRON_HEATING": 6, "POTENTIAL_MODEL": "WEIMER",
                "GPI_NCFILE": gpi, "KP": None, "JOULEFAC": 1.5, "ONEWAY": True,
                "GSWM_MI_DI_NCFILE": "/old/data/gswm_diurn_1.25d_99km.nc",
                "HE_COEFS_NCFILE": "/old/data/he_coefs_dres.nc", "OPDIFFCAP": "2e9",
                "other_input": [None]},
        "job": {"project_code": "P12345", "queue": "develop", "job_priority": "economy",
                "resource": {"select": 2, "ncpus": 128, "mpiprocs": 128}, "walltime": "02:00:00",
                "nprocs": 256, "other_job": ["module load foo"]},
    }


@pytest.fixture
def od():
    with open(tiegcmrun.OPTION_DESCRIPTIONS_FILE, encoding="utf-8") as f:
        return apply_system_config(json.load(f), socket.gethostname())


@pytest.fixture
def run(tmp_path, monkeypatch):
    """A saved run: <tmp>/run/stdout/tiegcmrun_parameters.json; TIEGCMDATA is an empty dir."""
    data = tmp_path / "data"
    data.mkdir()
    stdout = tmp_path / "run" / "stdout"
    stdout.mkdir(parents=True)
    monkeypatch.setenv("TIEGCMDATA", str(data))
    monkeypatch.setenv("TIEGCMHOME", str(tmp_path / "tiegcm"))
    path = stdout / tiegcmrun.PARAMETERS_JSON
    path.write_text(json.dumps(saved_parameters(str(data)), indent=4))
    return path


ASKED = {
    "simulation.job_name": "storm", "simulation.hpc_system": "derecho",
    "model.specification.horires": 1.25, "model.specification.zitop": 7,
    "inp.start_time": "2015-03-16T00:00:00", "inp.stop_time": "2015-03-19T00:00:00",
    "inp.solar_flux_level": "high", "inp.CALENDAR_ADVANCE": 1, "inp.NSTEP_SUB": 10,
    "inp.SECFLDS": ["TN", "UN", "VN"], "inp.ELECTRON_HEATING": 6, "inp.POTENTIAL_MODEL": "WEIMER",
    "inp.JOULEFAC": 1.5, "inp.ONEWAY": True,
    "job.derecho.project_code": "P12345", "job.derecho.queue": "develop",
    "job.derecho.job_priority": "economy", "job.derecho.walltime": "02:00:00",
    "job.derecho.other_job": ["module load foo"],
}


def test_load_seeds_asked_answers_only(run, od):
    seed, info = tiegcmrun.load_prefill(str(run), od)
    assert info == {"file": str(run), "ignored": [], "dropped": [], "notes": []}
    for key, value in ASKED.items():
        assert seed[key] == value, key
    # relative to the JSON's directory, as with -o
    assert seed["inp.GPI_NCFILE"] == str(run.parent.parent.parent / "data" / "gpi_test.nc")
    derived = ("inp.LABEL", "inp.PRIHIST", "inp.MXHIST_PRIM", "inp.SECHIST", "inp.OUTPUT",
               "inp.SECOUT", "inp.SOURCE", "inp.SOURCE_START", "inp.STEP", "inp.START_YEAR",
               "inp.secondary_start_time", "inp.GSWM_MI_DI_NCFILE", "inp.HE_COEFS_NCFILE",
               "inp.OPDIFFCAP", "model.specification.vertres", "model.specification.mres",
               "model.specification.nres_grid", "model.data.parentdir", "model.data.execdir", "model.data.workdir",
               "model.data.histdir", "model.data.modelexe", "model.data.log_file",
               "model.data.modeldir", "model.data.tgcmdata", "model.data.make",
               "model.data.input_file", "job.derecho.resource.select", "job.derecho.nprocs")
    assert not [k for k in derived if k in seed]
    # unset answers keep the default
    assert "inp.KP" not in seed and "inp.segment" not in seed and "inp.other_input" not in seed
    # only the file's machine block
    assert not [k for k in seed if k.startswith("job.") and not k.startswith("job.derecho.")]


def test_unknown_obsolete_and_invalid(run, od):
    saved = json.loads(run.read_text())
    saved["inp"]["NOT_A_KEY"] = 1
    saved["job"]["module_list"] = ["old"]
    saved["_comment"] = "kept quiet"
    saved["model"]["specification"]["horires"] = 3.0
    saved["simulation"]["job_name"] = "bad name"
    saved["inp"]["GPI_NCFILE"] = "/no/such/gpi.nc"
    run.write_text(json.dumps(saved))
    seed, info = tiegcmrun.load_prefill(str(run), od)
    assert info["ignored"] == ["inp.NOT_A_KEY", "job.module_list"]
    dropped = " ".join(info["dropped"])
    for key in ("model.specification.horires", "simulation.job_name", "inp.GPI_NCFILE"):
        assert key in dropped and key not in seed
    assert seed["inp.start_time"] == "2015-03-16T00:00:00"


def test_custom_deck_is_an_answer(run, od, tmp_path):
    deck = tmp_path / "my.inp"
    deck.write_text("&tgcm_input\n/\n")
    seed, _ = tiegcmrun.load_prefill({"model": {"data": {"input_file": str(deck)}}, "inp": None}, od)
    assert seed == {"model.data.input_file": str(deck)}


def test_job_keys_need_the_machine(od):
    seed, info = tiegcmrun.load_prefill({"job": {"queue": "main"}}, od)
    assert seed == {} and info["ignored"] == ["job.queue"]


def test_prefill_file_names(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert tiegcmrun.prefill_file("auto") is None
    (tmp_path / "tiegcmrun_input.json").write_text("{}")
    assert tiegcmrun.prefill_file("auto") == "tiegcmrun_input.json"
    assert tiegcmrun.prefill_file(str(tmp_path)) == str(tmp_path / "tiegcmrun_input.json")
    (tmp_path / "tiegcmrun_parameters.json").write_text("{}")
    assert tiegcmrun.prefill_file("auto") == "tiegcmrun_parameters.json"
    assert tiegcmrun.prefill_file(str(tmp_path)) == str(tmp_path / "tiegcmrun_parameters.json")
    with pytest.raises(ValueError):
        tiegcmrun.prefill_file(str(tmp_path / "missing.json"))
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="no tiegcmrun_parameters.json or tiegcmrun_input.json"):
        tiegcmrun.prefill_file(str(tmp_path / "empty"))


def test_template_input_json(od, tmp_path):
    t = tmp_path / "tiegcmrun_input.json"
    t.write_text(json.dumps({"inp": {"solar_flux_level": "low", "SECFLDS": ["TN", "NE"]}}))
    seed, info = tiegcmrun.load_prefill(str(tmp_path), od)
    assert seed == {"inp.solar_flux_level": "low", "inp.SECFLDS": ["TN", "NE"]}
    assert info["file"] == str(t)


def test_apply_prefill_defaults(od):
    seed = {"inp.segment": "1 0 0 0", "inp.JOULEFAC": 1.5, "job.derecho.queue": "develop"}
    before = copy.deepcopy(od)
    out = tiegcmrun.apply_prefill_defaults(od, seed)
    assert od == before                                  # a copy
    assert out["inp"]["segment"]["default"] == [1, 0, 0, 0]
    assert out["inp"]["JOULEFAC"]["default"] == 1.5 and out["inp"]["JOULEFAC"]["_prefilled"]
    assert out["job"]["derecho"]["queue"]["default"] == "develop"
    assert "_prefilled" not in out["inp"]["KP"]


def test_prompt_offers_and_keeps_the_saved_value(od, monkeypatch):
    out = tiegcmrun.apply_prefill_defaults(od, {"inp.solar_flux_level": "high", "inp.JOULEFAC": 1.5,
                                                "inp.segment": "1 0 0 0"})
    monkeypatch.setattr(builtins, "input", lambda prompt="": "")
    shown = []
    monkeypatch.setattr(builtins, "input", lambda prompt="": shown.append(prompt) or "")
    assert tiegcmrun.get_run_option("solar_flux_level", out["inp"]["solar_flux_level"]) == "high"
    assert "[\x1b[32mhigh\x1b[0m]" in shown[0]                # the GREEN default
    assert tiegcmrun.get_run_option("segment", out["inp"]["segment"]) == "1 0 0 0"
    # an EXPERT answer is used in BASIC too (not asked)
    assert tiegcmrun.get_run_option("JOULEFAC", out["inp"]["JOULEFAC"], "BASIC") == 1.5


def test_engage_and_build_pins_win():
    seed = {"inp.start_time": "x", "model.specification.horires": 1.25, "job.aitken.queue": "long",
            "inp.solar_flux_level": "high", "model.data.parentdir": "/old"}
    engage = {"start_time": "y", "horires": 2.5, "queue": "normal", "skip": ["parentdir"]}
    assert tiegcmrun.prefill_pinned(seed, engage) == {"inp.solar_flux_level": "high"}
    assert tiegcmrun.prefill_pinned(seed) == seed


def _args(**kw):
    base = dict(prefill=None, no_prefill=False, options_path=None, engage=None, onlycompile=False,
                verbose=False, benchmark=None, compile=False, mode="BASIC")
    base.update(kw)
    return argparse.Namespace(**base)


def test_auto_finds_nothing_prints_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    assert tiegcmrun._prefill_seed(_args()) == {}
    assert capsys.readouterr().out == ""


def test_auto_asks_first(run, monkeypatch, capsys):
    monkeypatch.chdir(run.parent)
    monkeypatch.setattr(tiegcmrun, "stdin_is_tty", lambda: True)    # pytest's stdin is not a terminal
    asked = []
    monkeypatch.setattr(builtins, "input", lambda prompt="": asked.append(prompt) or "n")
    assert tiegcmrun._prefill_seed(_args()) == {}
    assert asked == ["Found tiegcmrun_parameters.json: start the wizard from its values? [Y/n] "]
    assert capsys.readouterr().out == ""
    monkeypatch.setattr(builtins, "input", lambda prompt="": "")
    seed = tiegcmrun._prefill_seed(_args())
    assert seed["inp.solar_flux_level"] == "high"
    out = capsys.readouterr().out
    assert f"Prefilled {len(seed)} option(s) from tiegcmrun_parameters.json; derived values are recomputed." in out
    # --no-prefill never looks; -o never prefills
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    assert tiegcmrun._prefill_seed(_args(no_prefill=True)) == {}
    assert tiegcmrun._prefill_seed(_args(options_path="x.json")) == {}


def test_explicit_path_does_not_ask(run, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    seed = tiegcmrun._prefill_seed(_args(prefill=str(run), verbose=True))
    out = capsys.readouterr().out
    assert f"from {run};" in out and "  inp.solar_flux_level = high" in out
    assert seed["job.derecho.queue"] == "develop"


def test_cli_refuses_prefill_with_replay():
    with pytest.raises(ValueError, match="cannot be combined with -o"):
        tiegcmrun.tiegcmrun(["-o", "x.json", "--prefill", "y.json"])
    with pytest.raises(SystemExit):
        tiegcmrun.tiegcmrun(["--prefill", "--no-prefill"])


def _state(seed, **kw):
    tf = tiegcmrun._import_own_tui_form()
    return tf.build_state(_args(prefill_seed=seed, **kw), mode="BASIC")


def test_wizard_state_is_seeded(run, od):
    seed, _ = tiegcmrun.load_prefill(str(run), od)
    st = _state(seed)
    for key, value in ASKED.items():
        assert st.values[key] == value, key
        assert key in st.overridden
    # derived from the seeded answers, not taken from the file
    assert st.values["inp.LABEL"] == "storm_1.25x0.125"
    assert float(st.values["model.specification.mres"]) == 1.0
    # a new run directory: parentdir is not taken from the old run
    assert st.values["model.data.parentdir"] == "."
    assert st.values["model.data.execdir"] == os.path.join(".", "exec")
    assert st.values["inp.PRISTART"] == "75 0 0 0"
    assert "inp.PRIHIST" not in st.overridden
    # nothing seeded: the state is as before
    plain = _state({})
    assert not plain.overridden & set(seed)
    # a 'gen' answer is seeded as 'gen', not derived again
    assert _state({**seed, "inp.GPI_NCFILE": "gen"}).values["inp.GPI_NCFILE"] == "gen"


def test_wizard_engage_pins_win(run, od):
    seed, _ = tiegcmrun.load_prefill(str(run), od)
    engage = {"start_time": "2015-03-15T20:00:00", "stop_time": "2015-03-18T00:00:00",
              "segment": [0, 2, 0, 0], "horires": 2.5, "skip": [], "hpc_system": "derecho"}
    st = _state(seed, engage=engage)
    assert st.values["inp.start_time"] == "2015-03-15T20:00:00"
    assert st.values["model.specification.horires"] == 2.5
    assert st.values["inp.solar_flux_level"] == "high"


def test_build_state_accepts_a_dict():
    tf = tiegcmrun._import_own_tui_form()
    st = tf.build_state({"prefill_seed": {"inp.solar_flux_level": "medium"}}, mode="BASIC")
    assert st.values["inp.solar_flux_level"] == "medium"


def test_auto_not_a_tty(run, monkeypatch, capsys):
    monkeypatch.chdir(run.parent)
    monkeypatch.setattr(tiegcmrun, "stdin_is_tty", lambda: False)
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    assert tiegcmrun._prefill_seed(_args()) == {}
    out = capsys.readouterr().out
    assert "not a terminal" in out and "tiegcmrun_parameters.json" in out
    # a bare --prefill in a pipe is taken without asking
    assert tiegcmrun._prefill_seed(_args(prefill="auto"))["inp.solar_flux_level"] == "high"
    assert "Found" not in capsys.readouterr().out


def test_engage_does_not_ask_about_the_cwd_file(run, monkeypatch):
    monkeypatch.chdir(run.parent)
    monkeypatch.setattr(tiegcmrun, "stdin_is_tty", lambda: True)
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    assert tiegcmrun._prefill_seed(_args(engage={"skip": []})) == {}
    seed = tiegcmrun._prefill_seed(_args(prefill=str(run),
                                          engage={"skip": ["start_time"], "hpc_system": "derecho"}))
    assert seed["inp.solar_flux_level"] == "high" and "inp.start_time" not in seed


def test_machine_keyed_job_block(od):
    keyed = {"simulation": {"job_name": "x", "hpc_system": "aitken"},
             "job": {"_common": {"walltime": "01:00:00", "other_job": ["x"]},
                     "aitken": {"project_code": "s1234", "queue": "devel", "resource": {"model": "rom_ait"}},
                     "derecho": {"queue": "main"}}}
    flat = {"simulation": {"job_name": "x", "hpc_system": "aitken"},
            "job": {"walltime": "01:00:00", "other_job": ["x"], "project_code": "s1234", "queue": "devel",
                    "resource": {"model": "rom_ait"}}}
    seed, info = tiegcmrun.load_prefill(keyed, od)
    assert seed == tiegcmrun.load_prefill(flat, od)[0]
    assert seed["job.aitken.queue"] == "devel"
    assert info["ignored"] == []
    assert any("job.aitken is the job block" in n for n in info["notes"])
    assert tiegcmrun.load_prefill(keyed, od, hpc="derecho")[0]["job.derecho.queue"] == "main"


def test_bare_data_filename(run, od, tmp_path):
    data = tmp_path / "data"                       # the run fixture's $TIEGCMDATA
    gpi = data / "gpi_1960001-2025340.nc"
    gpi.write_text("")
    t = tmp_path / "tiegcmrun_input.json"
    t.write_text(json.dumps({"inp": {"GPI_NCFILE": "gpi_1960001-2025340.nc", "IMF_NCFILE": "imf_x.nc"}}))
    seed, info = tiegcmrun.load_prefill(str(t), od)
    assert seed["inp.GPI_NCFILE"] == str(gpi)
    assert "imf_x.nc" in " ".join(info["dropped"]) and "inp.IMF_NCFILE" not in seed
    (tmp_path / "imf_x.nc").write_text("")         # beside the JSON: relative to it, as with -o
    seed, info = tiegcmrun.load_prefill(str(t), od)
    assert seed["inp.IMF_NCFILE"] == str(tmp_path / "imf_x.nc")


def test_generated_forcing_seeds_gen(run, od):
    saved = json.loads(run.read_text())
    external = saved["inp"]["GPI_NCFILE"]
    gpi = run.parent / "gpi_2015075-2015078.nc"                       # in workdir
    (run.parent.parent / "hist").mkdir()
    imf = run.parent.parent / "hist" / "imf_2015075-2015078.nc"       # in parentdir
    gpi.write_text("")
    imf.write_text("")
    saved["inp"]["GPI_NCFILE"], saved["inp"]["IMF_NCFILE"] = str(gpi), str(imf)
    run.write_text(json.dumps(saved))
    seed, info = tiegcmrun.load_prefill(str(run), od)
    assert seed["inp.GPI_NCFILE"] == "gen" and seed["inp.IMF_NCFILE"] == "gen"
    assert len([n for n in info["notes"] if "was generated into the saved run" in n]) == 2
    # a file outside the saved run stays a path
    saved["inp"]["GPI_NCFILE"] = external
    run.write_text(json.dumps(saved))
    seed, _ = tiegcmrun.load_prefill(str(run), od)
    assert seed["inp.GPI_NCFILE"] == external


def test_external_executables_are_answers(run, od):
    saved = json.loads(run.read_text())
    saved["model"]["data"]["modelexe"] = "/builds/tiegcm_2.5/tiegcm.exe"
    saved["model"]["data"]["coupled_modelexe"] = "../exec/tiegcm.x"
    run.write_text(json.dumps(saved))
    seed, info = tiegcmrun.load_prefill(str(run), od)
    assert seed["model.data.modelexe"] == "/builds/tiegcm_2.5/tiegcm.exe"
    assert "model.data.coupled_modelexe" not in seed
    assert info["notes"] == []


def test_onlycompile_drops_the_executables(run, monkeypatch):
    saved = json.loads(run.read_text())
    saved["model"]["data"]["modelexe"] = "/builds/tiegcm_2.5/tiegcm.exe"
    run.write_text(json.dumps(saved))
    monkeypatch.setattr(builtins, "input", lambda prompt="": pytest.fail("asked"))
    assert tiegcmrun._prefill_seed(_args(prefill=str(run)))["model.data.modelexe"] == "/builds/tiegcm_2.5/tiegcm.exe"
    assert "model.data.modelexe" not in tiegcmrun._prefill_seed(_args(prefill=str(run), onlycompile=True))


def test_prefill_walltime_default():
    od = {"walltime": {"default": "24:00:00", "_prefilled": True}}
    note = tiegcmrun.prefill_walltime_default(od, "devel", "aitken")
    assert od["walltime"]["default"] == "02:00:00" and "_prefilled" not in od["walltime"]
    assert "more than queue 'devel'" in note
    od = {"walltime": {"default": "01:00:00", "_prefilled": True}}
    assert tiegcmrun.prefill_walltime_default(od, "devel", "aitken") is None
    assert od["walltime"]["default"] == "01:00:00"
    od = {"walltime": {"default": "08:00:00"}}                       # no --prefill: as upstream
    assert tiegcmrun.prefill_walltime_default(od, "devel", "aitken") is None
    assert od["walltime"]["default"] == "02:00:00"


def test_stale_forcing_default_dropped():
    desc = {"default": "/x/gpi_2015070-2015080.nc", "_prefilled": True}
    note = tiegcmrun.drop_stale_forcing_default(desc, "GPI_NCFILE", "2024-05-10T00:00:00", "2024-05-12T00:00:00")
    assert desc["default"] is None and "_prefilled" not in desc and "gpi_2015070-2015080.nc" in note
    desc = {"default": "/x/gpi_2015070-2015080.nc", "_prefilled": True}
    assert tiegcmrun.drop_stale_forcing_default(desc, "GPI_NCFILE", "2015-03-16T00:00:00",
                                                "2015-03-19T00:00:00") is None
    assert desc["default"] == "/x/gpi_2015070-2015080.nc"
    desc = {"default": "gen", "_prefilled": True}
    assert tiegcmrun.drop_stale_forcing_default(desc, "IMF_NCFILE", "2024-05-10T00:00:00",
                                                "2024-05-12T00:00:00") is None
    assert desc["default"] == "gen"


def test_forcing_cover_problem():
    import misc
    w = ("2024-05-10T00:00:00", "2024-05-12T00:00:00")
    assert misc.forcing_cover_problem("GPI_NCFILE", "gen", *w) is None
    assert misc.forcing_cover_problem("GPI_NCFILE", "/x/gpi_test.nc", *w) is None
    assert "enter 'gen'" in misc.forcing_cover_problem("GPI_NCFILE", "/x/gpi_2015070-2015080.nc", *w)
    assert misc.forcing_cover_problem("GPI_NCFILE", "/x/gpi_2024001-2024200.nc", *w) is None


def test_wizard_flags_a_stale_forcing_file(run, od, tmp_path):
    ext = tmp_path / "ext"
    ext.mkdir()
    gpi = ext / "gpi_2015070-2015080.nc"
    gpi.write_text("")
    saved = json.loads(run.read_text())
    saved["inp"]["GPI_NCFILE"] = str(gpi)
    run.write_text(json.dumps(saved))
    seed, _ = tiegcmrun.load_prefill(str(run), od)
    assert seed["inp.GPI_NCFILE"] == str(gpi)
    st = _state(seed)
    assert "inp.GPI_NCFILE" not in st.errors()               # the saved window is covered
    st.set("inp.start_time", "2024-05-10T00:00:00")
    st.set("inp.stop_time", "2024-05-12T00:00:00")
    assert "covers 2015070..2015080" in st.errors()["inp.GPI_NCFILE"]
    st.set("inp.GPI_NCFILE", "gen")
    assert "inp.GPI_NCFILE" not in st.errors()
    # any dated answer is checked, not only a prefilled one
    plain = _state({})
    plain.set("inp.start_time", "2024-05-10T00:00:00")
    plain.set("inp.stop_time", "2024-05-12T00:00:00")
    plain.set("inp.GPI_NCFILE", str(gpi))
    assert "covers 2015070..2015080" in plain.errors()["inp.GPI_NCFILE"]


def test_reprompt_while_re_asks_until_covered(monkeypatch, capsys):
    desc = {"LEVEL": "BASIC", "prompt": "GPI file", "default": None}
    asked = []
    def fake(name, description, mode="BASIC", skip_parameters=[], run_dir=None):
        asked.append((name, run_dir))
        return "gen"
    monkeypatch.setattr(tiegcmrun, "get_run_option", fake)
    import misc
    w = ("2024-05-10T00:00:00", "2024-05-12T00:00:00")
    out = tiegcmrun._reprompt_while("GPI_NCFILE", desc, "BASIC", [], "/x/gpi_2015070-2015080.nc",
                                    lambda v: misc.forcing_cover_problem("GPI_NCFILE", v, *w), run_dir="/r")
    assert out == "gen" and asked == [("GPI_NCFILE", "/r")]
    assert "covers 2015070..2015080" in capsys.readouterr().out
    # not prompted at this mode: refused, not re-asked
    with pytest.raises(ValueError, match="covers 2015070..2015080"):
        tiegcmrun._reprompt_while("GPI_NCFILE", desc, "BASIC", ["GPI_NCFILE"],
                                  "/x/gpi_2015070-2015080.nc",
                                  lambda v: misc.forcing_cover_problem("GPI_NCFILE", v, *w))
    assert len(asked) == 1
