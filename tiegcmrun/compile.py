import os
import shutil
import subprocess
import filecmp
import sys
from textwrap import dedent
import logging

def compile_tiegcm(options, debug, coupling = False, hidra = False):
    """
    Compiles the TIEGCM model with the given options.

    Args:
        options (dict): A dictionary containing the model options.
        debug (bool): A boolean indicating whether to enable debug mode.
        coupling (bool): A boolean indicating whether to enable coupling.

    Returns:None
    """
    o = options
    modeldir  = o["model"]["data"]["modeldir"]
    execdir   = o["model"]["data"]["execdir"]
    workdir = o["model"]["data"]["workdir"]
    outdir = o["model"]["data"]["histdir"]
    tgcmdata  = o["model"]["data"]["tgcmdata"]
    utildir   = os.path.join(o["model"]["data"]["modeldir"],"scripts")
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
    nres_grid = float(o["model"]["specification"]["nres_grid"])
    make      = o["model"]["data"]["make"]
    coupling  = coupling
    hidra     = hidra

    if coupling == True:
        modelexe = os.path.basename(o["model"]["data"]["coupled_modelexe"])
        model = o["model"]["data"]["coupled_modelexe"]
    else:
        modelexe = os.path.basename(o["model"]["data"]["modelexe"])
        model = o["model"]["data"]["modelexe"]
    debug = debug
    
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

    # Set srcdir based on modeldir
    srcdir = os.path.join(modeldir, 'src')

    # Check if srcdir exists
    if not os.path.isdir(srcdir):
        print(f">>> Cannot find model source directory {srcdir} <<<")
        sys.exit(1)

    # Convert srcdir to an absolute path
    srcdir = os.path.abspath(srcdir)  

    if tgcmdata == "None":
        tgcmdata = os.environ['TIEGCMDATA']
        print(f"Set tgcmdata = {tgcmdata}")

    if not os.path.isdir(tgcmdata):
        print(f">>> Cannot find data directory {tgcmdata}")

    # Check horizontal resolution
    if horires not in [5, 2.5, 1.25, 0.625]:
        print(f">>> Unknown model horizontal resolution {horires} <<<")
        sys.exit(1)
    
    # Check vertical resolution
    if vertres not in [0.5, 0.25, 0.125, 0.0625]:
        print(f">>> Unknown model vertical resolution {vertres} <<<")
        sys.exit(1)
    
    if nres_grid == "None" or nres_grid == None:
        if mres == 2:
            nres_grid = 5
        elif mres == 1:
            nres_grid = 6
        elif mres == 0.5:
            nres_grid = 7
        else:
            print(f">>> Unsupported magnetic resolution {mres} <<<")
            sys.exit(1)

    # Copy make if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, os.path.basename(make))):
        shutil.copy(os.path.join(utildir, os.path.basename(make)), execdir)
    # Copy Makefile if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, 'Makefile')):
        shutil.copy(os.path.join(utildir, 'Makefile'), execdir)
    # Copy mkdepends if it does not exist in execdir
    if not os.path.isfile(os.path.join(execdir, 'mkdepends')):
        shutil.copy(os.path.join(utildir, 'mkdepends'), execdir)
    
    if input == '' or output == '':
        input = os.path.abspath(input)
        output = os.path.abspath(output)
        
    util = os.path.abspath(utildir)


    coupling_file_path = os.path.join(execdir, 'coupling')

    # Check if the coupling file exists
    if os.path.isfile(coupling_file_path):
        with open(coupling_file_path, 'r') as file:
            lastcoupling = file.read().strip().lower() == 'true'
        # Compare coupling values
        if lastcoupling != coupling:
            print(f"Clean execdir {execdir} because coupling flag switched from {lastcoupling} to {coupling}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)
            with open(coupling_file_path, 'w') as file:
                file.write(str(coupling))
    else:
        # Create the coupling file and write the coupling value
        with open(coupling_file_path, 'w') as file:
            file.write(str(coupling))
        print(f"Created file coupling with coupling flag = {coupling}")


    hidra_file_path = os.path.join(execdir, 'hidra')

    # Check if the hidra file exists
    if os.path.isfile(hidra_file_path):
        with open(hidra_file_path, 'r') as file:
            lasthidra = file.read().strip().lower() == 'true'
        # Compare hidra values
        if lasthidra != hidra:
            print(f"Clean execdir {execdir} because hidra flag switched from {lasthidra} to {hidra}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)
            with open(hidra_file_path, 'w') as file:
                file.write(str(hidra))
    else:
        # Create the hidra file and write the hidra value
        with open(hidra_file_path, 'w') as file:
            file.write(str(hidra))
        print(f"Created file hidra with hidra flag = {hidra}")

    debug_file_path = os.path.join(execdir, 'debug')

    # Check if the debug file exists
    if os.path.isfile(debug_file_path):
        with open(debug_file_path, 'r') as file:
            lastdebug = file.read().strip().lower() == 'true'

        # Compare debug values
        if lastdebug != debug:
            print(f"Clean execdir {execdir} because debug flag switched from {lastdebug} to {debug}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)

            with open(debug_file_path, 'w') as file:
                file.write(str(debug))
    else:
        # Create the debug file and write the debug value
        with open(debug_file_path, 'w') as file:
            file.write(str(debug))
        print(f"Created file debug with debug flag = {debug}")


    # Create the defs.h content
    defs_content = dedent(f"""\
    #define DLAT {horires}
    #define DLON {horires}
    #define GLON1 -180
    #define DLEV {vertres}
    #define ZIBOT -7
    #define ZITOP {zitop}
    #define NRES_GRID {nres_grid}
    """)

    # Write to defs.h
    defs_path = 'defs.h'
    with open(defs_path, 'w') as file:
        file.write(defs_content)

    # Check if defs.h exists in execdir and compare
    execdir_defs_path = os.path.join(execdir, 'defs.h')
    if os.path.isfile(execdir_defs_path):
        if not filecmp.cmp(defs_path, execdir_defs_path, shallow=False):
            # Files differ, switch resolutions
            print(f"Switching defs.h for model resolution {horires} x {vertres}")
            mycwd = os.getcwd()
            os.chdir(execdir)
            subprocess.run(['gmake', 'clean'])
            os.chdir(mycwd)
            shutil.copy(defs_path, execdir_defs_path)
        else:
            print(f"defs.h already set for model resolution {horires} x {vertres}")
    else:
        # defs.h does not exist in execdir, copy the file
        print(f"Copying {defs_path} to {execdir_defs_path} for resolution {horires} x {vertres}")
        shutil.copy(defs_path, execdir_defs_path)


    try:
        os.chdir(execdir)
        print(f"\nBegin building {model} in {os.getcwd()}")
    except OSError:
        print(f">>> Cannot cd to execdir {execdir}")
        sys.exit(1)



    # Create Make.env file
    make_env_path = os.path.join(execdir, 'Make.env')
    with open(make_env_path, 'w') as file:
        file.write(f"""MAKE_MACHINE  = {make}
DIRS          = . {srcdir}
EXECNAME      = {model}
NAMELIST      = {input}
OUTPUT        = {output}
COUPLING      = {str(coupling).upper()}
HIDRA         = {str(hidra).upper()}
DEBUG         = {str(debug).upper()}
""")

    # Build the model
    try:
        subprocess.run(['gmake', '-j8', 'all'], check=True)
        shutil.copy(model, workdir)
        print(f"Executable copied from {model} to {workdir}")
    except subprocess.CalledProcessError:
        print(">>> Error return from gmake all")
        sys.exit(1)