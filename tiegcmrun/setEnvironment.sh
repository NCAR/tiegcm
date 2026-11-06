#!/usr/bin/env bash
# usage: source setEnvironment.sh [--machine <name>]
# The machine (config/machines.yaml nodename_regex matched against the host name, or --machine /
# SYSTEM), its modules and its conda activation come from tiegcmrun/config/machines.yaml: the same
# entry the job scripts use. A host that matches no cluster loads no modules.
# TIEGCM_SETENV_DRYRUN=1: print the module and conda commands instead of running them.
TIEGCMHOME=""
TIEGCMDATA=""
SYSTEM="" # a config/machines.yaml machine name (empty: from the host name)
SAVED_MODULES="" # Name of existing set of modules saved on system
CONDA_ENV="" # Name of existing conda env

if [ "${1:-}" = "--machine" ] && [ -n "${2:-}" ]; then
    SYSTEM="$2"
fi

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-${(%):-%x}}" )" &> /dev/null && pwd )"
TIEGCM_MACHINE_ENV="$SCRIPT_DIR/machine_env.py"

# Run one shell command (print it under TIEGCM_SETENV_DRYRUN=1).
_tiegcm_run() {
    if [ -n "${TIEGCM_SETENV_DRYRUN:-}" ]; then
        printf '+ %s\n' "$1"
    else
        eval "$1"
    fi
}

# Run each line of machine_env.py's output for this machine.
_tiegcm_run_machine() {
    local out line
    if ! out="$(python3 "$TIEGCM_MACHINE_ENV" --machine "$SYSTEM" "$@")"; then
        echo "WARNING: could not read $SCRIPT_DIR/config/machines.yaml (python3 with PyYAML needed): skipped '$*'"
        return 1
    fi
    TIEGCM_RAN=0
    while IFS= read -r line; do
        [ -n "$line" ] && _tiegcm_run "$line" && TIEGCM_RAN=$((TIEGCM_RAN + 1))
    done <<< "$out"
    return 0
}

if [ -z "$SYSTEM" ]; then
    SYSTEM="$(python3 "$TIEGCM_MACHINE_ENV" --hostname "${HOSTNAME:-$(hostname)}" name)"
fi
if [ -z "$SYSTEM" ]; then
    echo "WARNING: no machine from $SCRIPT_DIR/config/machines.yaml (python3 with PyYAML needed): no modules loaded"
else
    echo "System set: $SYSTEM"
fi

if [ -z "$TIEGCMHOME" ]; then
    ROOT_DIR="$(echo "$SCRIPT_DIR" | sed 's:/tiegcmrun$::')"
    TIEGCMHOME="$ROOT_DIR"
    read -p "Enter tiegcm model directory [$TIEGCMHOME]:" TIEGCMHOME_INP
    if [ -z "$TIEGCMHOME_INP" ]; then
        echo "TIEGCMHOME set: $TIEGCMHOME"
        export TIEGCMHOME="$TIEGCMHOME"
    else
        echo "TIEGCMHOME set: $TIEGCMHOME_INP"
        export TIEGCMHOME="$TIEGCMHOME_INP"
    fi
else
    echo "TIEGCMHOME set: $TIEGCMHOME"
fi


if [ -z "$TIEGCMDATA" ]; then
    read -p "Enter tiegcm data directory [$TIEGCMDATA]:" TIEGCMDATA_INP
    export TIEGCMDATA="$TIEGCMDATA_INP"
    echo "TIEGCMDATA set: $TIEGCMDATA_INP"
else
    export TIEGCMDATA="$TIEGCMDATA"
    echo "TIEGCMDATA set: $TIEGCMDATA"
fi

if [ -z "$SAVED_MODULES" ]; then
    if [ -n "$SYSTEM" ] && _tiegcm_run_machine modules && [ "$TIEGCM_RAN" -gt 0 ]; then
        _tiegcm_run "module list"
        read -p "Save Module List (y/n)? " answer
        if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
            read -p "Enter name of module list: " SAVED_MODULESINP
            SAVED_MODULES="$SAVED_MODULESINP"
        fi
    fi
else
    _tiegcm_run "module restore $SAVED_MODULES"
    echo "Restore Modules: $SAVED_MODULES"
    _tiegcm_run "module list"
fi

REQUIREMENTSTXT="${TIEGCMHOME}/tiegcmrun/requirements.txt"

# conda lines that need no environment name (e.g. module load conda)
[ -n "$SYSTEM" ] && _tiegcm_run_machine conda-setup

_tiegcm_activate() {
    if [ -n "$SYSTEM" ]; then
        _tiegcm_run_machine conda "$1"
    else
        _tiegcm_run "conda activate \"$1\""
    fi
}

if [ -z "$CONDA_ENV" ]; then
    read -p "Enter conda environment name: " CONDA_ENVINP
    if [ -n "${TIEGCM_SETENV_DRYRUN:-}" ] && [ -n "$CONDA_ENVINP" ]; then
        echo "Environment '$CONDA_ENVINP' (not checked: dry run)."
    elif [ -n "$CONDA_ENVINP" ] && conda env list | grep -q "^${CONDA_ENVINP} "; then
        echo "Environment '$CONDA_ENVINP' found."
    else
        echo "Environment '$CONDA_ENVINP' does not exist."
        CONDA_ENVINP=""
    fi

    if [ -z "$CONDA_ENVINP" ]; then
        read -p "Enter new conda environment name: " CONDA_ENVINP
        PKGS_DIRS=""
        [ -n "$SYSTEM" ] && PKGS_DIRS="$(python3 "$TIEGCM_MACHINE_ENV" --machine "$SYSTEM" conda-pkgs-dirs)"
        if [ -n "$PKGS_DIRS" ]; then
            _tiegcm_run "export CONDA_PKGS_DIRS=$PKGS_DIRS"
        fi
        _tiegcm_run "conda create -n \"$CONDA_ENVINP\" python=3.8"
    fi
else
    CONDA_ENVINP="$CONDA_ENV"
fi
_tiegcm_activate "$CONDA_ENVINP"
_tiegcm_run "pip install -r \"${REQUIREMENTSTXT}\""
