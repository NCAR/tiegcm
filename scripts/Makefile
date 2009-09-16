#
#-----------------------------------------------------------------------
#
# Makefile for building TIEGCM
# Supported platforms include IBM-AIX, SGI-IRIX, Linux
#
#-----------------------------------------------------------------------
#
# The following variables are added to this template with the appropriate
# values by tgcm_config. If tgcm_config is not called, or the Makefile 
# needs to be debugged or customized, the variables can be set manually here.
#
# UNAMES     := [operating system from "uname -s" command]
# UNAMEN     := [machine name from "uname -n" command]
# EXECNAME   := [executable name (usually model name with version number)]
# MODEL_CLASS:= [model class name (tiegcm, timegcm, glbmean, mtgcm, etc)]
# EXECDIR    := [executable directory]
# UTILDIR    := [location of script utilities, Makefile, tgcm_config, etc]
# LIB_NETCDF := [location of netcdf library libnetcdf.a]
# INC_NETCDF := [location of netcdf fortran include file]
# DEBUG      := [TRUE or FALSE]
# MPI        := [TRUE if running parallel code]
# MSS        := [TRUE if NCAR mass store is available]
# NBIT       := [32 or 64 bit Linux machine]
#
null  :=
space := $(null) $(null)
comma := $(null),$(null)
#
#-----------------------------------------------------------------------
# Load dependency search paths.
#
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH    := $(subst $(space),:,$(VPATH))               
#------------------------------------------------------------------------
# Default goal (force Mkdepends before making exec):
#
all: Depends 
#
Depends::
	perl $(UTILDIR)/Mkdepends
	$(MAKE) -f $(MAKENAME) $(EXECDIR)/$(EXECNAME)

CURDIR   := $(shell pwd)
MAKENAME := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
SOURCES  := $(shell cat Srcfiles)
OBJS     := $(addsuffix .o, $(basename $(SOURCES)))
LIBS     := -L$(LIB_NETCDF) -lnetcdf
LDFLAGS  :=
#
# This is tiegcm stand-alone, so no coupling with cism:
COUPLING := NO 
#
# Note FFLAGS and CPPFLAGS are included by gmake.
#
$(EXECDIR)/$(EXECNAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)
#-----------------------------------------------------------------------
#
# IBM operating system:
#
ifeq ($(UNAMES),AIX)
  CPPFLAGS := -WF -DAIX
  ifeq ($(MPI),TRUE)
    CPPFLAGS += -DMPI
  endif
  ifeq ($(MSS),TRUE)
    CPPFLAGS += -DMSS
  endif
  ifeq ($(COUPLING), INTERCOMM)
    CPPFLAGS += -DINTERCOMM
  endif
  ifeq ($(COUPLING), ADHOC)
    CPPFLAGS += -DCISMAH
  endif
  OPTIM    := -O3 -qstrict
  ifeq ($(DEBUG),TRUE)
    CPPFLAGS += -DDEBUG
    OPTIM    := -g
    DBGFLAGS := -qinitauto=7FF7FFFF -qflttrap=zero:inv:en -qsigtrap=xl__trcedump
#   DBGFLAGS := -C -qzerosize -qinitauto=7FF7FFFF -qflttrap=zero:inv:en \
#               -qsigtrap=xl__trcedump
  else
    LIBS     += -lmass       # IBM optimized lib
  endif
  ifeq ($(MPI),TRUE)
     FC         := mpxlf_r
   else
     FC         := xlf90
  endif
  CPPFLAGS := $(subst $(space),$(comma),$(CPPFLAGS))
  FFLAGS     := $(cpp_path) $(OPTIM) $(DBGFLAGS) -qarch=auto -qrealsize=8 -q64
#
# 11/21/08: Add 64K page sizes:
  LDFLAGS    += -bloadmap:loadmap -q64 $(OPTIM)
  LDFLAGS    += -bdatapsize:64K -bstackpsize:64K -btextpsize:64K
endif # AIX
#-----------------------------------------------------------------------
#
# SGI operating system:
#
ifeq ($(UNAMES),IRIX64)
  CPPFLAGS := -DIRIX
  ifeq ($(MPI),TRUE)
    CPPFLAGS += -DMPI
  endif
  ifeq ($(MSS),TRUE)
    CPPFLAGS += -DMSS
  endif
  OPTIM    := -O3
  ifeq ($(DEBUG),TRUE)
    CPPFLAGS += -DDEBUG
    OPTIM    := -g
  endif
  FC         := f90
  FFLAGS     := $(cpp_path) $(OPTIM) -r8 -macro_expand -OPT:Olimit=4479 -v 
  LIBS       += -L/usr/lib32/mips4 -lffio -lcomplib.sgimath -lblas_mp -lfpe
  ifeq ($(MPI),TRUE)
    LIBS       += -L/lusr/local/lib32/r4i4 -lmpi
  endif
  LDFLAGS    += -mips4 -r10000
endif
#-----------------------------------------------------------------------
#
# Linux operating system:
#
ifeq ($(UNAMES),Linux)
#
# Sets for all linux systems:
#
  CPPFLAGS := -DLINUX
  ifeq ($(MPI),TRUE)
    CPPFLAGS += -DMPI
  endif
  ifeq ($(MSS),TRUE)
    CPPFLAGS += -DMSS
  endif
  ifeq ($(COUPLING), INTERCOMM)
    CPPFLAGS += -DINTERCOMM
  endif
  ifeq ($(COUPLING), ADHOC)
    CPPFLAGS += -DCISMAH
  endif
  OPTIM    := -O2
  ifeq ($(DEBUG),TRUE)
    CPPFLAGS += -DDEBUG
    OPTIM    := -g
  endif
#
# Machine specific sets, within Linux OS:
#
# Default is pgf90 (e.g., serial linux systems at HAO):
#
  FC       := pgf90
  CPPFLAGS += -I/usr/local/include
  FFLAGS   := $(cpp_path) $(OPTIM) -r8
  ifeq ($(MPI),TRUE)
    LIBS   += -L/usr/local/lib -lmpich
  endif
#
# 64-bit netcdf4 hdf libs:
  ifeq ($(NBIT),64)
    LIBS += -L$(LIB_NETCDF) -lhdf5_hl -lhdf5 -lsz -L/usr/lib64 -lz
#
# Default is pgf90 (e.g., serial linux systems at HAO),
# or mpif90 if MPI is set:
#
    ifeq ($(MPI),TRUE)
      FC := mpif90
    endif
  endif
#
# SCD Linux system "gale" (serial):
#
  ifeq ($(UNAMEN),gale)                     # CISL/DASG Linux "gale" (no mpi)
    FC       := ifort
    FFLAGS := $(cpp_path) $(OPTIM) -r8 -fixed -diag-disable remark
  endif
#
# SCD Linux cluster "lightning" (with or without MPI):
#
  ifeq ($(UNAMEN),lightning)               # CISL Linux cluster "lightning"
    ifeq ($(MPI),TRUE)
      FC := mpif90
    else
      FC := pgf90
    endif
    cpp_path += -I/usr/local/mpich/include
    FFLAGS := $(cpp_path) $(OPTIM) -r8
    LIBS := -L$(LIB_NETCDF) -lnetcdf -L/usr/lib64 -lz
    LIBS += -L/usr/local/mpich/lib -lmpich
  endif
endif # linux
#-----------------------------------------------------------------------
#
RM := rm
clean:
	$(RM) -f Depends Srcfiles *.o *.mod $(EXECDIR)/$(EXECNAME)
#
# Include object file targets and dependencies made by Mkdepends:
#
include Depends

