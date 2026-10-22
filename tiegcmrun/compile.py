import os
import shutil
import subprocess
import filecmp
import sys
from textwrap import dedent
from misc import mres_to_nres_grid, default_make_file, git_state, svn_revision_stamp, tiegcm_env

def _gmake_clean(execdir):
    """`gmake clean` in execdir; exits on failure so stale objects are never linked."""
    rc = subprocess.run(['gmake', 'clean'], cwd=execdir).returncode
    if rc != 0:
        print(f">>> Error return {rc} from gmake clean in {execdir}; not building on stale objects <<<")
        sys.exit(1)


def compile_tiegcm(options, debug, coupling = False, hidra = False):
    """Build TIE-GCM, restoring the caller's working directory afterwards."""
    cwd = os.getcwd()
    try:
        _compile_tiegcm(options, debug, coupling, hidra)
    finally:
        os.chdir(cwd)


def _compile_tiegcm(options, debug, coupling = False, hidra = False):
    """Compile TIE-GCM; `debug` builds a Fortran debug executable."""
    o = options
    # absolute, since the build chdirs
    modeldir  = os.path.abspath(o["model"]["data"]["modeldir"])
    execdir   = os.path.abspath(o["model"]["data"]["execdir"])
    workdir = os.path.abspath(o["model"]["data"]["workdir"])
    outdir = os.path.abspath(o["model"]["data"]["histdir"])
    tgcmdata  = o["model"]["data"]["tgcmdata"]
    utildir   = os.path.join(modeldir,"scripts")
    try:
        input     = o["model"]["data"]["input_file"]
    except:
        input = ""
    try:
        output    = o["model"]["data"]["log_file"]
    except:
        output = ""
    horires   = float(o["model"]["specification"]["horires"])
    vertres   = float(o["model"]["specification"]["vertres"])
    zitop     = float(o["model"]["specification"]["zitop"])
    mres      = float(o["model"]["specification"]["mres"])
    nres_grid = o["model"]["specification"].get("nres_grid")
    nres_grid = float(nres_grid) if nres_grid is not None else float(mres_to_nres_grid(mres))
    make      = o["model"]["data"].get("make") or default_make_file(
        modeldir, o["simulation"]["hpc_system"])
    if not make:
        print(f">>> No Make fragment for hpc_system {o['simulation']['hpc_system']}: "
              f"set model.data.make <<<")
        sys.exit(1)
    coupling  = coupling
    hidra     = hidra

    if coupling == True:
        modelexe = os.path.basename(o["model"]["data"]["coupled_modelexe"])
        model = os.path.abspath(o["model"]["data"]["coupled_modelexe"])
    else:
        modelexe = os.path.basename(o["model"]["data"]["modelexe"])
        model = os.path.abspath(o["model"]["data"]["modelexe"])
    debug = debug
    input = os.path.abspath(input) if input else input
    output = os.path.abspath(output) if output else output

    try:
        os.makedirs(workdir)
    except:
        print(f"{workdir} exitsts")
    try:
        os.makedirs(outdir)
    except:
        print(f"{outdir} exitsts")
    try:
        os.makedirs(execdir)
    except:
        print(f"{execdir} exitsts")
    os.chdir(workdir)

    if not os.path.isdir(modeldir):
        print(f">>> Cannot find model directory {modeldir} <<<")
        sys.exit(1)

    if not os.path.isdir(utildir):
        print(f">>> Cannot find model directory {utildir} <<<")
        sys.exit(1)

    srcdir = os.path.join(modeldir, 'src')

    if not os.path.isdir(srcdir):
        print(f">>> Cannot find model source directory {srcdir} <<<")
        sys.exit(1)

    srcdir = os.path.abspath(srcdir)  

    if tgcmdata == "None":
        tgcmdata = tiegcm_env("TIEGCMDATA")
        print(f"Set tgcmdata = {tgcmdata}")

    if not os.path.isdir(tgcmdata):
        print(f">>> Cannot find data directory {tgcmdata}")

    if horires not in [5, 2.5, 1.25, 0.625]:
        print(f">>> Unknown model horizontal resolution {horires} <<<")
        sys.exit(1)

    if vertres not in [0.5, 0.25, 0.125, 0.0625]:
        print(f">>> Unknown model vertical resolution {vertres} <<<")
        sys.exit(1)

    # refresh execdir's copy of the Make fragment whenever the source changed
    make_src = make if os.path.isfile(make) else os.path.join(utildir, os.path.basename(make))
    make_copy = os.path.join(execdir, os.path.basename(make))
    if os.path.isfile(make_src):
        if not os.path.isfile(make_copy) or not filecmp.cmp(make_src, make_copy, shallow=False):
            shutil.copy(make_src, make_copy)
    elif not os.path.isfile(make_copy):
        print(f">>> Cannot find Make fragment {make} <<<")
        sys.exit(1)
    if not os.path.isfile(os.path.join(execdir, 'Makefile')):
        shutil.copy(os.path.join(utildir, 'Makefile'), execdir)
    if not os.path.isfile(os.path.join(execdir, 'mkdepends')):
        shutil.copy(os.path.join(utildir, 'mkdepends'), execdir)
    
    util = os.path.abspath(utildir)


    # flag files in execdir record the last build's settings; a change forces gmake clean
    coupling_file_path = os.path.join(execdir, 'coupling')

    if os.path.isfile(coupling_file_path):
        with open(coupling_file_path, 'r') as file:
            lastcoupling = file.read().strip().lower() == 'true'
        if lastcoupling != coupling:
            print(f"Clean execdir {execdir} because coupling flag switched from {lastcoupling} to {coupling}")
            _gmake_clean(execdir)
            with open(coupling_file_path, 'w') as file:
                file.write(str(coupling))
    else:
        with open(coupling_file_path, 'w') as file:
            file.write(str(coupling))
        print(f"Created file coupling with coupling flag = {coupling}")


    hidra_file_path = os.path.join(execdir, 'hidra')

    if os.path.isfile(hidra_file_path):
        with open(hidra_file_path, 'r') as file:
            lasthidra = file.read().strip().lower() == 'true'
        if lasthidra != hidra:
            print(f"Clean execdir {execdir} because hidra flag switched from {lasthidra} to {hidra}")
            _gmake_clean(execdir)
            with open(hidra_file_path, 'w') as file:
                file.write(str(hidra))
    else:
        with open(hidra_file_path, 'w') as file:
            file.write(str(hidra))
        print(f"Created file hidra with hidra flag = {hidra}")

    debug_file_path = os.path.join(execdir, 'debug')

    if os.path.isfile(debug_file_path):
        with open(debug_file_path, 'r') as file:
            lastdebug = file.read().strip().lower() == 'true'

        if lastdebug != debug:
            print(f"Clean execdir {execdir} because debug flag switched from {lastdebug} to {debug}")
            _gmake_clean(execdir)

            with open(debug_file_path, 'w') as file:
                file.write(str(debug))
    else:
        with open(debug_file_path, 'w') as file:
            file.write(str(debug))
        print(f"Created file debug with debug flag = {debug}")


    defs_content = dedent(f"""\
    #define DLAT {horires}
    #define DLON {horires}
    #define GLON1 -180
    #define DLEV {vertres}
    #define ZIBOT -7
    #define ZITOP {zitop}
    #define NRES_GRID {nres_grid}
    """)

    defs_path = 'defs.h'
    with open(defs_path, 'w') as file:
        file.write(defs_content)

    execdir_defs_path = os.path.join(execdir, 'defs.h')
    if os.path.isfile(execdir_defs_path):
        if not filecmp.cmp(defs_path, execdir_defs_path, shallow=False):
            print(f"Switching defs.h for model resolution {horires} x {vertres}")
            _gmake_clean(execdir)
            shutil.copy(defs_path, execdir_defs_path)
        else:
            print(f"defs.h already set for model resolution {horires} x {vertres}")
    else:
        print(f"Copying {defs_path} to {execdir_defs_path} for resolution {horires} x {vertres}")
        shutil.copy(defs_path, execdir_defs_path)


    try:
        os.chdir(execdir)
        print(f"\nBegin building {model} in {os.getcwd()}")
    except OSError:
        print(f">>> Cannot cd to execdir {execdir}")
        sys.exit(1)



    # MAKE_MACHINE is the bare name of the execdir copy: `make veryclean` deletes $(MAKE_MACHINE).
    # SVN_REVISION is the model's git short hash, stored in every history (16 chars).
    tiegcm_git = git_state(modeldir)
    make_env_path = os.path.join(execdir, 'Make.env')
    with open(make_env_path, 'w') as file:
        file.write(f"""MAKE_MACHINE  = {os.path.basename(make)}
DIRS          = . {srcdir}
EXECNAME      = {model}
NAMELIST      = {input}
OUTPUT        = {output}
COUPLING      = {str(coupling).upper()}
HIDRA         = {str(hidra).upper()}
DEBUG         = {str(debug).upper()}
SVN_REVISION  = {svn_revision_stamp(tiegcm_git)}
""")

    try:
        subprocess.run(['gmake', '-j8', 'all'], check=True)
    except subprocess.CalledProcessError:
        print(">>> Error return from gmake all")
        sys.exit(1)
    dest = os.path.join(workdir, os.path.basename(model))
    if os.path.exists(dest) and os.path.samefile(model, dest):
        print(f"Executable {model} is already in {workdir}")
    else:
        shutil.copy(model, workdir)
        print(f"Executable copied from {model} to {workdir}")
