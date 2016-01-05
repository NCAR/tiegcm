# ESMF application makefile fragment
#
# Use the following ESMF_ variables to compile and link
# your ESMF application against this ESMF build.
#
# !!! VERY IMPORTANT: If the location of this ESMF build is   !!!
# !!! changed, e.g. libesmf.a is copied to another directory, !!!
# !!! this file - esmf.mk - must be edited to adjust to the   !!!
# !!! correct new path                                        !!!
#
# Please see end of file for options used on this ESMF build
#


ESMF_VERSION_STRING=6.3.0r

ESMF_VERSION_MAJOR=6
ESMF_VERSION_MINOR=3
ESMF_VERSION_REVISION=0
ESMF_VERSION_PATCHLEVEL=0
ESMF_VERSION_PUBLIC='T'
ESMF_VERSION_BETASNAPSHOT='F'


ESMF_APPSDIR=/home/tgcm/esmf/pgi-9.04/apps/appsO/Linux.pgi.64.mpich.default
ESMF_LIBSDIR=/home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default


ESMF_F90COMPILER=/opt/local/pgi-9.04/linux86-64/9.0/mpi/mpich/bin/mpif90
ESMF_F90LINKER=/opt/local/pgi-9.04/linux86-64/9.0/mpi/mpich/bin/mpif90

ESMF_F90COMPILEOPTS=-O -fpic -mcmodel=small -lpthread
ESMF_F90COMPILEPATHS=-I/home/tgcm/esmf/pgi-9.04/mod/modO/Linux.pgi.64.mpich.default -I/home/tgcm/esmf/pgi-9.04/src/include 
ESMF_F90COMPILECPPFLAGS=-DESMF_NO_INTEGER_1_BYTE -DESMF_NO_INTEGER_2_BYTE -DESMF_LAPACK=1 -DESMF_MOAB=1 -DESMF_MPIIO -DESMF_NO_OPENMP -DESMF_NO_OPENACC -DSx86_64_small=1 -DESMF_OS_Linux=1 -DESMF_MPICH -DESMF_PGIVERSION_MAJOR=9 -DESMF_PGIVERSION_MINOR=0 -DESMF_PGIVERSION_PATCH=4
ESMF_F90COMPILEFREECPP=
ESMF_F90COMPILEFREENOCPP=-Mfreeform
ESMF_F90COMPILEFIXCPP=-Mpreprocess -Mnofreeform
ESMF_F90COMPILEFIXNOCPP=

ESMF_F90LINKOPTS= -mcmodel=small -lpthread -Wl,--no-as-needed
ESMF_F90LINKPATHS=-L/home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default  -L/opt/local/pgi-9.04/linux86-64/9.0-4/lib
ESMF_F90LINKRPATHS=-Wl,-rpath,/home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default  -Wl,-rpath,/opt/local/pgi-9.04/linux86-64/9.0-4/lib
ESMF_F90LINKLIBS= -lpmpich++ -lmpich -pgcpplibs -ldl
ESMF_F90ESMFLINKLIBS=-lesmf  -lpmpich++ -lmpich -pgcpplibs -ldl

ESMF_CXXCOMPILER=mpiCC
ESMF_CXXLINKER=mpiCC

ESMF_CXXCOMPILEOPTS=-O -DNDEBUG -fpic -mcmodel=small -lpthread
ESMF_CXXCOMPILEPATHS= -I/home/tgcm/esmf/pgi-9.04/src/include  
ESMF_CXXCOMPILECPPFLAGS=-DESMF_NO_INTEGER_1_BYTE -DESMF_NO_INTEGER_2_BYTE -DESMF_LAPACK=1 -DESMF_MOAB=1 -DESMF_MPIIO -DESMF_NO_OPENMP -DESMF_NO_OPENACC -DSx86_64_small=1 -DESMF_OS_Linux=1 -D__SDIR__='' -DESMF_MPICH

ESMF_CXXLINKOPTS= -mcmodel=small -lpthread -Wl,--no-as-needed
ESMF_CXXLINKPATHS=-L/home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default  -L/opt/local/pgi-9.04/linux86-64/9.0-4/lib
ESMF_CXXLINKRPATHS=-Wl,-rpath,/home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default  -Wl,-rpath,/opt/local/pgi-9.04/linux86-64/9.0-4/lib
ESMF_CXXLINKLIBS= -lmpich -pgf90libs -ldl
ESMF_CXXESMFLINKLIBS=-lesmf  -lmpich -pgf90libs -ldl

ESMF_SO_F90COMPILEOPTS=-fpic
ESMF_SO_F90LINKOPTS=-shared
ESMF_SO_F90LINKOPTSEXE=-Wl,-export-dynamic
ESMF_SO_CXXCOMPILEOPTS=-fpic
ESMF_SO_CXXLINKOPTS=-shared
ESMF_SO_CXXLINKOPTSEXE=-Wl,-export-dynamic

ESMF_OPENMP_F90COMPILEOPTS= -mp
ESMF_OPENMP_F90LINKOPTS= -mp
ESMF_OPENMP_CXXCOMPILEOPTS= -mp --exceptions
ESMF_OPENMP_CXXLINKOPTS= -mp --exceptions

ESMF_OPENACC_F90COMPILEOPTS= -acc -Minfo
ESMF_OPENACC_F90LINKOPTS= -acc -Minfo
ESMF_OPENACC_CXXCOMPILEOPTS= -acc -Minfo
ESMF_OPENACC_CXXLINKOPTS= -acc -Minfo

# Internal ESMF variables, do NOT depend on these!

ESMF_INTERNAL_DIR=/home/tgcm/esmf/pgi-9.04

#
# !!! The following options were used on this ESMF build !!!
#
# ESMF_DIR: /home/tgcm/esmf/pgi-9.04
# ESMF_OS: Linux
# ESMF_MACHINE: x86_64
# ESMF_ABI: 64
# ESMF_COMPILER: pgi
# ESMF_BOPT: O
# ESMF_COMM: mpich
# ESMF_SITE: default
# ESMF_PTHREADS: ON
# ESMF_OPENMP: OFF
# ESMF_OPENACC: OFF
# ESMF_ARRAY_LITE: FALSE
# ESMF_NO_INTEGER_1_BYTE: TRUE
# ESMF_NO_INTEGER_2_BYTE: TRUE
# ESMF_FORTRANSYMBOLS: default
# ESMF_DEFER_LIB_BUILD:   ON
# 
# ESMF environment variables pointing to 3rd party software:
# ESMF_MOAB:              internal
# ESMF_LAPACK:            internal
