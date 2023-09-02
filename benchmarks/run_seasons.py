
import os
import subprocess

# Setting up variables
home = "/glade/u/home/nikhilr/tiegcm_func"
modelroot = f"{home}/tiegcm2.0"
tgcmdata = "/glade/p/hao/tgcm/data/tiegcm2.0"
tgcmrun = f"{modelroot}/tgcmrun/tgcmrun"
queue = "regular"
submit = "yes"
execute = "yes"
resolutions = ['5.0', '2.5']
runs = [
    'decsol_smax', 'decsol_smin', 'junsol_smax', 'junsol_smin',
    'mareqx_smax', 'mareqx_smin', 'sepeqx_smax', 'sepeqx_smin'
]

resolutions = ['2.5']
runs = [
    'decsol_smin']

# Looping through resolutions and runs
for res in resolutions:
    case = f"tiegcm_res{res}_seasons"
    work = f"{home}/{case}"
    
    for run in runs:
        print(f"------------------ {run} --------------------")
        stdout = f"{work}/tiegcm_res{res}_{run}/stdout"
        histdir = f"{work}/tiegcm_res{res}_{run}/hist"
        execdir = f"{work}/tiegcm_res{res}_{run}/exec"
        
        # Command execution
        cmd = [
            tgcmrun,
            '-q', queue,
            '-run_name', run,
            '-model_name', 'tiegcm',
            '-model_res', res,
            '-model_root', modelroot,
            '-stdout_dir', stdout,
            '-execdir', execdir,
            '-hist_dir', histdir,
            '-submit', submit,
            '-tgcmdata', tgcmdata,
            '-execute', execute
        ]
        subprocess.run(cmd)

