#!/usr/bin/env bash
# usage: setupEnvironment.sh

TIEGCMHOME="/home7/nrao3/tiegcm"
TIEGCMDATA="/nobackup/nrao3/tiegcm/data/tiegcm3.0"
SYSTEM="" # derecho or pleiades
SAVED_MODULES="tiegcm" # Name of existing set of modules saved on system 
CONDA_ENV="my_tiegcm" # Name of existing conda env


if [ -z "$SYSTEM" ]; then
    if [[ $HOSTNAME == *"derecho"* ]]; then
        SYSTEM="derecho"
        echo "System set: $SYSTEM"
    elif [[ $HOSTNAME == *"pfe"* ]]; then
        SYSTEM="pleiades"
        echo "System set: $SYSTEM"
    else
        SYSTEM="custom"
    fi
else
    echo "System set: $SYSTEM"
fi

if [ -z "$TIEGCMHOME" ]; then
    SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-${(%):-%x}}" )" &> /dev/null && pwd )"
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
    if [ "$SYSTEM" == "derecho" ]; then
        module purge
        module load conda
        module load ncarenv/23.06
        module load intel/2023.0.0
        module load cray-mpich/8.1.25
        module load mkl/2023.0.0
        module load netcdf-mpi/4.9.2
        module load craype/2.7.20
        module load ncarcompilers/1.0.0
        module load hdf5-mpi/1.12.2
        module load esmf/8.5.0
        
        module list

        read -p "Save Module List (y/n)? " answer
        if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
            read -p "Enter name of module list:" SAVED_MODULESINP
            SAVED_MODULES="$SAVED_MODULESINP"
        fi
    elif [ "$SYSTEM" == "pleiades" ]; then
        module purge
        module load nas
        module load pkgsrc/2023Q3
        module load comp-intel/2020.4.304 
        module load mpi-hpe/mpt.2.28_25Apr23_rhel87
        module load mpi-hpe/mpt 
        module load szip/2.1.1
        module load hdf4/4.2.12
        module load hdf5/1.8.18_mpt
        module load netcdf/4.4.1.1_mpt
        module use -a /swbuild/analytix/tools/modulefiles
        module load miniconda3/v4
        module list
        
        read -p "Save Module List (y/n)? " answer
        if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
            read -p "Enter name of module list: " SAVED_MODULESINP
            SAVED_MODULES="$SAVED_MODULESINP"
        fi
    fi
else
    module restore $SAVED_MODULES
    echo "Restore Modules: $SAVED_MODULES"
    module list
fi

REQUIREMENTSTXT="${TIEGCMHOME}/tiegcmrun/requirements.txt"

if [ -z "$CONDA_ENV" ]; then
    read -p "Enter conda environment name: " CONDA_ENVINP
    if conda env list | grep -q "^${CONDA_ENVINP} "; then
        echo "Environment '$CONDA_ENVINP' found."
    else
        echo "Environment '$CONDA_ENVINP' does not exist."
        CONDA_ENVINP=""
    fi

    if [ -z "$CONDA_ENVINP" ]; then
        read -p "Enter new conda environment name: " CONDA_ENVINP
        conda create -n "$CONDA_ENVINP" python=3.8
        if [ "$SYSTEM" == "pleiades" ]; then
            export CONDA_PKGS_DIRS=/nobackup/$USER/.conda/pkgs
        fi
        conda activate "$CONDA_ENVINP"
        pip install -r "${REQUIREMENTSTXT}"
    else
        if [ "$SYSTEM" == "pleiades" ]; then
            source activate "$CONDA_ENVINP"
        else
            conda activate "$CONDA_ENVINP"

        fi
        pip install -r "${REQUIREMENTSTXT}"
        # INSTALLREQ=FALSE
        # trap '' PIPE
        # while IFS= read -r package || [[ -n "$package" ]]; do
        #    if ! pip list | grep -q "^${package}[[:space:]]"; then
        #        echo "Package '$package' is not installed in the '$CONDA_ENVINP' environment."
        #        INSTALLREQ=TRUE
        #        break
        #    fi
        # done < "$REQUIREMENTSTXT"
        # if [ "$INSTALLREQ" == "TRUE" ]; then
        #    pip install -r "${REQUIREMENTSTXT}"
        # fi
    fi
else
    if [ "$SYSTEM" == "pleiades" ]; then
        source activate "$CONDA_ENV"
    else
        conda activate "$CONDA_ENV"
    fi
    pip install -r "${REQUIREMENTSTXT}"
fi



