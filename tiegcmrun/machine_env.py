#!/usr/bin/env python3
"""The shell environment of a config/machines.yaml machine: module and conda lines.

Used by setEnvironment.sh (with the system python3: needs only PyYAML) and by the job scripts
(output_solver). The machine is the one whose nodename_regex matches the host name, or --machine.

  machine_env.py [--machine M] [--hostname H] name             the machine's name
  machine_env.py [--machine M] [--hostname H] modules          module purge/use/load lines (+ local_prefix, bash)
  machine_env.py [--machine M] [--hostname H] conda-setup      conda lines that need no environment name
  machine_env.py [--machine M] [--hostname H] conda ENV        the lines that activate conda environment ENV
                                                               (a venv when ENV/bin/activate exists)
  machine_env.py [--machine M] [--hostname H] conda-pkgs-dirs  CONDA_PKGS_DIRS for conda create ('' if none)
"""
import argparse
import copy
import os
import socket
import sys

from jobgen import JobGen

CONFIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")


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


def machine_module_lines(jg):
    """module purge/use/load lines of the machine's own modules: a stacks block in tiegcmrun's
    machines.yaml holds coupled-run facts (stacks.grt.status), never the tiegcm job's modules."""
    if not jg.m.get("stacks"):
        return jg.module_lines()
    own = copy.copy(jg)
    own.m = {k: v for k, v in jg.m.items() if k != "stacks"}
    return own.module_lines()


def shell_module_lines(jg):
    """The module lines of an interactive bash session (build and tools), local_prefix included."""
    lines = machine_module_lines(jg)
    if jg.m.get("local_prefix"):
        lines += local_prefix_lines(jg.m["local_prefix"], "bash")
    return lines


def conda_lines(jg):
    """The machine's bash conda activation lines (conda_activate_sh, else conda_activate)."""
    return list(jg.m.get("conda_activate_sh") or jg.m.get("conda_activate") or [])


def venv_lines(jg):
    """The machine's bash venv activation lines (venv_activate_sh, else venv_activate)."""
    return list(jg.m.get("venv_activate_sh") or jg.m.get("venv_activate") or [])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--machine", default=None, help="a config/machines.yaml machine name")
    p.add_argument("--hostname", default=None, help="the host name to resolve (default: this host)")
    p.add_argument("what", choices=("name", "modules", "conda-setup", "conda", "conda-pkgs-dirs"))
    p.add_argument("env", nargs="?", default=None, help="conda: the environment name, or a venv path (a directory holding bin/activate)")
    a = p.parse_args(argv)
    try:
        jg = JobGen(CONFIG_DIR, machine=a.machine or None, hostname=a.hostname or socket.gethostname())
    except ValueError as e:
        print(f"machine_env.py: {e}", file=sys.stderr)
        return 2
    if a.what == "name":
        out = [jg.name]
    elif a.what == "modules":
        out = shell_module_lines(jg)
    elif a.what == "conda-setup":
        out = [ln for ln in conda_lines(jg) if "{conda_env}" not in ln]
    elif a.what == "conda":
        if not a.env:
            p.error("conda needs the environment name")
        venv = os.path.isfile(os.path.join(a.env, "bin", "activate"))
        if venv and not venv_lines(jg):
            print(f"machine_env.py: config/machines.yaml has no venv_activate lines for machine "
                  f"{jg.name!r}: cannot activate {a.env}.", file=sys.stderr)
            return 2
        out = [ln.replace("{conda_env}", a.env) for ln in (venv_lines(jg) if venv else conda_lines(jg))
               if "{conda_env}" in ln]
    else:
        out = [str(jg.m.get("conda_pkgs_dirs") or "")]
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
