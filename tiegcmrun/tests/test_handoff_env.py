"""The hand-off job activates coupling.conda_env: a venv path (<env>/bin/activate exists) through
machines.yaml venv_activate, a conda name through conda_activate."""
import os
import shutil
import sys

import pytest
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

import machine_env  # noqa: E402
import output_solver  # noqa: E402
from output_solver import handoff_job_lines  # noqa: E402


def make_venv(tmp_path):
    venv = tmp_path / "myvenv"
    (venv / "bin").mkdir(parents=True)
    (venv / "bin" / "activate").write_text("# bash activate\n")
    (venv / "bin" / "activate.csh").write_text("# csh activate\n")
    return str(venv)


def config_without(tmp_path, machine, key):
    """A copy of config/ whose machine entry lacks key."""
    cfg = tmp_path / "config"
    shutil.copytree(output_solver.CONFIG_DIR, cfg)
    with open(cfg / "machines.yaml") as f:
        data = yaml.safe_load(f)
    data["machines"][machine].pop(key, None)
    data["machines"][machine].pop(key + "_sh", None)
    with open(cfg / "machines.yaml", "w") as f:
        yaml.safe_dump(data, f)
    return str(cfg)


def test_conda_name_uses_conda_activate():
    assert handoff_job_lines("aitken", "kaipy", "/w/handoff.py") == [
        "source /swbuild/analytix/tools/miniconda3_220407/bin/activate.csh kaipy",
        "python /w/handoff.py"]
    assert handoff_job_lines("derecho", "kaipy", "/w/handoff.py") == [
        "module load conda", "conda activate kaipy", "python /w/handoff.py"]


def test_venv_path_uses_venv_activate_csh(tmp_path):
    venv = make_venv(tmp_path)
    for machine in ("aitken", "derecho"):           # tcsh job shells
        assert handoff_job_lines(machine, venv, "/w/handoff.py") == [
            'if (! $?prompt) set prompt = ""', f"source {venv}/bin/activate.csh", "python /w/handoff.py"]


def test_venv_path_uses_venv_activate_bash(tmp_path):
    venv = make_venv(tmp_path)
    for machine in ("linux", "perlmutter"):         # bash job shells
        assert handoff_job_lines(machine, venv, "/w/handoff.py") == [
            f"source {venv}/bin/activate", "python /w/handoff.py"]


def test_missing_directory_is_a_conda_name(tmp_path):
    # a path without bin/activate is not a venv: the conda lines are kept (conda accepts a prefix)
    path = str(tmp_path / "no-such-env")
    assert handoff_job_lines("derecho", path, "s.py") == [
        "module load conda", f"conda activate {path}", "python s.py"]


def test_machine_without_venv_activate_refuses(tmp_path, monkeypatch):
    venv = make_venv(tmp_path)
    monkeypatch.setattr(output_solver, "CONFIG_DIR", config_without(tmp_path, "aitken", "venv_activate"))
    with pytest.raises(ValueError) as e:
        handoff_job_lines("aitken", venv, "/w/handoff.py")
    assert str(e.value) == ("config/machines.yaml has no venv_activate lines for machine 'aitken': "
                            "cannot write the job lines that run /w/handoff.py.")


def test_machine_without_conda_activate_refuses(tmp_path, monkeypatch):
    monkeypatch.setattr(output_solver, "CONFIG_DIR", config_without(tmp_path, "aitken", "conda_activate"))
    with pytest.raises(ValueError) as e:
        handoff_job_lines("aitken", "kaipy", "/w/handoff.py")
    assert str(e.value) == ("config/machines.yaml has no conda_activate lines for machine 'aitken': "
                            "cannot write the job lines that run /w/handoff.py.")


def test_every_machine_has_venv_activate():
    with open(os.path.join(output_solver.CONFIG_DIR, "machines.yaml")) as f:
        machines = yaml.safe_load(f)["machines"]
    for name, m in machines.items():
        assert m.get("venv_activate"), name
        csh = m.get("shell") in ("tcsh", "csh")
        assert m["venv_activate"][-1] == "source {conda_env}/bin/activate" + (".csh" if csh else ""), name


def test_machine_env_conda_venv_bash_lines(tmp_path, capsys):
    venv = make_venv(tmp_path)
    assert machine_env.main(["--machine", "aitken", "conda", venv]) == 0
    assert capsys.readouterr().out == f"source {venv}/bin/activate\n"
    assert machine_env.main(["--machine", "aitken", "conda", "kaipy"]) == 0
    assert capsys.readouterr().out == "source activate kaipy\n"


def test_machine_env_conda_venv_missing_lines(tmp_path, monkeypatch, capsys):
    venv = make_venv(tmp_path)
    monkeypatch.setattr(machine_env, "CONFIG_DIR", config_without(tmp_path, "aitken", "venv_activate"))
    assert machine_env.main(["--machine", "aitken", "conda", venv]) == 2
    assert capsys.readouterr().err == (f"machine_env.py: config/machines.yaml has no venv_activate lines "
                                       f"for machine 'aitken': cannot activate {venv}.\n")
