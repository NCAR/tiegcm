#!/bin/bash

# This script automates the build process for ESMF, netcdf-c, and netcdf-fortran on the Atiken HPC system.

purge_modules() {
  echo -e "\n\n\n\n\n\n\n\n\nPurging all loaded modules..."
  module --force purge
}



load_modules() {
  echo -e "\n\n\n\n\n\n\n\n\nConfiguring module paths and loading required modules..."
  module use -a /nasa/modulefiles/testing
  module use -a /swbuild/analytix/tools/modulefiles
  module load comp-intel/2020.4.304
  module load mpi-hpe/mpt.2.30
  module load hdf5/$H5_version
  module load miniconda3/v4
}

set_compiler_env() {
  echo -e "\n\n\n\n\n\n\n\n\nSetting compiler environment variables..."
  export FC=mpif90
  export CC=mpicc
  export CXX=mpicxx
}

set_dirs() {
  echo -e "\n\n\n\n\n\n\n\n\nSetting PREFIX and WORK directories..."
  export PREFIX=$HOME/local_atiken
  export WORK=$HOME/build_atiken
}

ensure_dirs() {
  echo -e "\n\n\n\n\n\n\n\n\nEnsuring PREFIX and WORK directories exist..."
  if [ ! -d "$PREFIX" ]; then
    mkdir -p "$PREFIX"
  fi
  if [ ! -d "$WORK" ]; then
    mkdir -p "$WORK"
  fi
}



download_tarballs() {
  echo -e "\n\n\n\n\n\n\n\n\nChecking for ESMF tarball..."
  if [ ! -f v${esmf_version}.tar.gz ]; then
    wget https://github.com/esmf-org/esmf/archive/refs/tags/v${esmf_version}.tar.gz
  fi
  echo -e "\n\n\n\n\n\n\n\n\nChecking for netcdf-c tarball..."
  if [ ! -f netcdf-c-${netcdf_c_version}.tar.gz ]; then
    wget https://downloads.unidata.ucar.edu/netcdf-c/${netcdf_c_version}/netcdf-c-${netcdf_c_version}.tar.gz
  fi
  echo -e "\n\n\n\n\n\n\n\n\nChecking for netcdf-fortran tarball..."
  if [ ! -f netcdf-fortran-${netcdf_fortran_version}.tar.gz ]; then
    wget https://downloads.unidata.ucar.edu/netcdf-fortran/${netcdf_fortran_version}/netcdf-fortran-${netcdf_fortran_version}.tar.gz
  fi
}

extract_tarballs() {
  echo -e "\n\n\n\n\n\n\n\n\nExtracting ESMF..."
  if [ -d "esmf-${esmf_version}" ]; then
    rm -rf esmf-${esmf_version}
  fi
  tar -xvzf v${esmf_version}.tar.gz
  echo -e "\n\n\n\n\n\n\n\n\nExtracting netcdf-c..."
  if [ -d "netcdf-c-${netcdf_c_version}" ]; then
    rm -rf netcdf-c-${netcdf_c_version}
  fi
  tar -xvzf netcdf-c-${netcdf_c_version}.tar.gz
  echo -e "\n\n\n\n\n\n\n\n\nExtracting netcdf-fortran..."
  if [ -d "netcdf-fortran-${netcdf_fortran_version}" ]; then
    rm -rf netcdf-fortran-${netcdf_fortran_version}
  fi
  tar -xvzf netcdf-fortran-${netcdf_fortran_version}.tar.gz
}

set_build_env() {
  export H5DIR=/nasa/hdf5/$H5_version
  export CFLAGS="-O2 -fPIC"
  export CXXFLAGS="-O2 -fPIC"
  export FCFLAGS="-O2 -fPIC"
  export CPPFLAGS="-I$H5DIR/include"
  export LDFLAGS="-L$H5DIR/lib"
  export LD_LIBRARY_PATH="$PREFIX/lib:$H5DIR/lib:$LD_LIBRARY_PATH"
  export PATH="$PREFIX/bin:$PATH"
}

build_netcdf_c() {
  echo -e "\n\n\n\n\n\n\n\n\nBuilding netcdf-c..."
  cd $WORK/netcdf-c-${netcdf_c_version}
  ./configure --prefix="$PREFIX" --with-hdf5="$H5DIR" --enable-netcdf-4 --enable-shared --enable-parallel-tests \
    --disable-dap \
    --disable-byterange \
    --disable-libxml2 \
    --disable-nczarr-zip
  echo -e "\n\n\n\n\n\n\n\n\nCompiling netcdf-c..."
  make -j 8
  echo -e "\n\n\n\n\n\n\n\n\nInstalling netcdf-c..."
  make install
  cd ..
}

build_netcdf_fortran() {
  echo -e "\n\n\n\n\n\n\n\n\nBuilding netcdf-fortran..."
  cd $WORK/netcdf-fortran-${netcdf_fortran_version}
  export CC=mpicc FC=mpif90 CPPFLAGS="-I$PREFIX/include $CPPFLAGS" LDFLAGS="-L$PREFIX/lib $LDFLAGS"
  echo -e "\n\n\n\n\n\n\n\n\nConfiguring netcdf-fortran..."
  ./configure --prefix="$PREFIX" --enable-shared
  echo -e "\n\n\n\n\n\n\n\n\nCompiling netcdf-fortran..."
  make -j 8
  echo -e "\n\n\n\n\n\n\n\n\nInstalling netcdf-fortran..."
  make install
  cd ..
}

build_esmf() {
  echo -e "\n\n\n\n\n\n\n\n\nBuilding ESMF..."
  cd $WORK/esmf-${esmf_version}
  export ESMF_DIR=$(pwd)
  export ESMF_COMPILER=intelgcc
  export ESMF_COMM=mpt
  #export ESMF_NETCDF=nc-config
  echo -e "\n\n\n\n\n\n\n\n\nCompiling ESMF..."
  make
  echo -e "\n\n\n\n\n\n\n\n\nInstalling ESMF..."
  export ESMF_INSTALL_PREFIX=$PREFIX/esmf
  make install
  cd ..
}

# Main script execution

purge_modules

echo -e "\n\n\n\n\n\n\n\n\nSetting environment variables for versions..."
export H5_version="1.12.3_mpt"
export netcdf_c_version="4.9.3"
export netcdf_fortran_version="4.6.2"
export esmf_version="8.9.0"

load_modules

set_compiler_env

set_dirs
ensure_dirs

echo -e "\n\n\n\n\n\n\n\n\nChanging to WORK directory: $WORK"
cd $WORK

download_tarballs
extract_tarballs
set_build_env
build_netcdf_c
build_netcdf_fortran
build_esmf




# export FC=mpif90
# export CC=mpicc
# export CXX=mpicxx

# export PREFIX=$HOME/local_atiken_test
# export LIBRARY_PATH=${LIBRARY_PATH}:$PREFIX/lib
# export LD_LIBRARY_PATH=$LIBRARY_PATH
# export CPATH=$PREFIX/include
# export PATH=${PATH}:$PREFIX/bin