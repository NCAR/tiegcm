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


ESMF_VERSION_STRING=6.3.0rp1

ESMF_VERSION_MAJOR=6
ESMF_VERSION_MINOR=3
ESMF_VERSION_REVISION=0
ESMF_VERSION_PATCHLEVEL=1
ESMF_VERSION_PUBLIC='T'
ESMF_VERSION_BETASNAPSHOT='F'


ESMF_APPSDIR=/home/tgcm/esmf/intel/esmf_6_3_0rp1/apps/appsO/Linux.intel.64.intelmpi.default
ESMF_LIBSDIR=/home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default


ESMF_F90COMPILER=mpiifort
ESMF_F90LINKER=mpiifort

ESMF_F90COMPILEOPTS=-O -fPIC -m64 -mcmodel=small -threads  -openmp
ESMF_F90COMPILEPATHS=-I/home/tgcm/esmf/intel/esmf_6_3_0rp1/mod/modO/Linux.intel.64.intelmpi.default -I/home/tgcm/esmf/intel/esmf_6_3_0rp1/src/include 
ESMF_F90COMPILECPPFLAGS=-DESMF_NO_INTEGER_1_BYTE -DESMF_NO_INTEGER_2_BYTE -DESMF_LAPACK=1 -DESMF_MOAB=1 -DESMF_MPIIO -DESMF_NO_OPENACC -DSx86_64_small=1 -DESMF_OS_Linux=1
ESMF_F90COMPILEFREECPP=
ESMF_F90COMPILEFREENOCPP=
ESMF_F90COMPILEFIXCPP=
ESMF_F90COMPILEFIXNOCPP=

ESMF_F90LINKOPTS= -m64 -mcmodel=small -threads -Wl,--no-as-needed  -openmp
ESMF_F90LINKPATHS=-L/home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default  -L/opt/local/intel2011/impi/4.0.1.007/intel64/lib -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64 -L/opt/local/intel2011/composerxe-2011.0.084/mkl/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -L/lib64 -L/lib -L/usr/lib64 -L/usr/lib
ESMF_F90LINKRPATHS=-Wl,-rpath,/home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default 
ESMF_F90LINKLIBS= -ldl -ldl -ldl -ldl -lpthread -lpthread -lpthread -lpthread -lrt -limf -lsvml -lm -lipgo -ldecimal -liomp5 -lcilkrts -lstdc++ -lgcc_s -lgcc -lirc -lpthread -lgcc_s -lgcc -lirc_s -ldl -lrt -ldl
ESMF_F90ESMFLINKLIBS=-lesmf  -ldl -ldl -ldl -ldl -lpthread -lpthread -lpthread -lpthread -lrt -limf -lsvml -lm -lipgo -ldecimal -liomp5 -lcilkrts -lstdc++ -lgcc_s -lgcc -lirc -lpthread -lgcc_s -lgcc -lirc_s -ldl -lrt -ldl

ESMF_CXXCOMPILER=mpiicpc
ESMF_CXXLINKER=mpiicpc

ESMF_CXXCOMPILEOPTS=-O -DNDEBUG -fPIC -m64 -mcmodel=small -pthread  -openmp
ESMF_CXXCOMPILEPATHS= -I/home/tgcm/esmf/intel/esmf_6_3_0rp1/src/include  
ESMF_CXXCOMPILECPPFLAGS=-DESMF_NO_INTEGER_1_BYTE -DESMF_NO_INTEGER_2_BYTE -DESMF_LAPACK=1 -DESMF_MOAB=1 -DESMF_MPIIO -DESMF_NO_OPENACC -DSx86_64_small=1 -DESMF_OS_Linux=1 -D__SDIR__=''

ESMF_CXXLINKOPTS= -m64 -mcmodel=small -pthread -Wl,--no-as-needed  -openmp
ESMF_CXXLINKPATHS=-L/home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default  -L/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64/
ESMF_CXXLINKRPATHS=-Wl,-rpath,/home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default  -Wl,-rpath,/opt/local/intel2011/composerxe-2011.0.084/compiler/lib/intel64/
ESMF_CXXLINKLIBS= -ldl -ldl -ldl -ldl -lpthread -lpthread -lpthread -lpthread -lrt -lifport -lifcoremt -limf -lsvml -lm -lipgo -liomp5 -lirc -lpthread -lgcc_s -lgcc -lirc_s -ldl -lrt -ldl
ESMF_CXXESMFLINKLIBS=-lesmf  -ldl -ldl -ldl -ldl -lpthread -lpthread -lpthread -lpthread -lrt -lifport -lifcoremt -limf -lsvml -lm -lipgo -liomp5 -lirc -lpthread -lgcc_s -lgcc -lirc_s -ldl -lrt -ldl

ESMF_SO_F90COMPILEOPTS=-fPIC
ESMF_SO_F90LINKOPTS=-shared
ESMF_SO_F90LINKOPTSEXE=-Wl,-export-dynamic
ESMF_SO_CXXCOMPILEOPTS=-fPIC
ESMF_SO_CXXLINKOPTS=-shared
ESMF_SO_CXXLINKOPTSEXE=-Wl,-export-dynamic

ESMF_OPENMP_F90COMPILEOPTS= -openmp
ESMF_OPENMP_F90LINKOPTS= -openmp
ESMF_OPENMP_CXXCOMPILEOPTS= -openmp
ESMF_OPENMP_CXXLINKOPTS= -openmp

ESMF_OPENACC_F90COMPILEOPTS=
ESMF_OPENACC_F90LINKOPTS=
ESMF_OPENACC_CXXCOMPILEOPTS=
ESMF_OPENACC_CXXLINKOPTS=

# Internal ESMF variables, do NOT depend on these!

ESMF_INTERNAL_DIR=/home/tgcm/esmf/intel/esmf_6_3_0rp1

#
# !!! The following options were used on this ESMF build !!!
#
# ESMF_DIR: /home/tgcm/esmf/intel/esmf_6_3_0rp1
# ESMF_OS: Linux
# ESMF_MACHINE: x86_64
# ESMF_ABI: 64
# ESMF_COMPILER: intel
# ESMF_BOPT: O
# ESMF_COMM: intelmpi
# ESMF_SITE: default
# ESMF_PTHREADS: ON
# ESMF_OPENMP: ON
# ESMF_OPENACC: OFF
# ESMF_ARRAY_LITE: FALSE
# ESMF_NO_INTEGER_1_BYTE: TRUE
# ESMF_NO_INTEGER_2_BYTE: TRUE
# ESMF_FORTRANSYMBOLS: default
# ESMF_DEFER_LIB_BUILD: ON
# ESMF_SHARED_LIB_BUILD: ON
# 
# ESMF environment variables pointing to 3rd party software:
# ESMF_MOAB:              internal
# ESMF_LAPACK:            internal
