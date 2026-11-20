"""A run generated on this host for another cluster: data files resolved here are named in a WARNING."""
import os
import socket
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

from replay_solver import host_resolved_files_warning  # noqa: E402


def options_for(hpc_system, data_dir):
    return {"simulation": {"hpc_system": hpc_system},
            "model": {"data": {"tgcmdata": "/glade/x"}, "specification": {"horires": 2.5}},
            "inp": {"SOURCE": f"{data_dir}/junsol_f70.nc", "GPI_NCFILE": f"{data_dir}/gpi.nc",
                    "HE_COEFS_NCFILE": None}}


def test_host_resolved_files_warning_cross_host(monkeypatch, tmp_path):
    monkeypatch.setattr(socket, "gethostname", lambda: "pfe27")
    monkeypatch.setenv("TIEGCMDATA", str(tmp_path))
    text = host_resolved_files_warning(options_for("derecho", tmp_path))
    assert text is not None
    assert text.startswith("WARNING: simulation.hpc_system is derecho but this host is aitken: ")
    assert f"inp.SOURCE, inp.GPI_NCFILE resolved under this host's TIEGCMDATA ({tmp_path})" in text
    assert "inp.HE_COEFS_NCFILE not found under derecho's tgcmdata (/glade/x)" in text
    assert text.endswith("; the derecho jobs will not find them: set them to derecho paths in the "
                         "tiegcmrun options JSON, or generate the run on derecho.")


def test_host_resolved_files_warning_cross_host_paths_elsewhere(monkeypatch, tmp_path):
    monkeypatch.setattr(socket, "gethostname", lambda: "pfe27")
    monkeypatch.setenv("TIEGCMDATA", str(tmp_path))
    options = options_for("derecho", "/glade/x")
    options["inp"]["HE_COEFS_NCFILE"] = "/glade/x/he.nc"
    assert host_resolved_files_warning(options) is None


def test_host_resolved_files_warning_same_host(monkeypatch, tmp_path):
    monkeypatch.setattr(socket, "gethostname", lambda: "pfe27")
    monkeypatch.setenv("TIEGCMDATA", str(tmp_path))
    assert host_resolved_files_warning(options_for("aitken", tmp_path)) is None
    assert host_resolved_files_warning(options_for("linux", tmp_path)) is None
