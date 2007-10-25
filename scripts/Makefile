#
#-----------------------------------------------------------------------
#
# Makefile for building TGCM (tiegcm or timegcm) 
# Supported platforms include IBM-AIX, SGI-IRIX, Linux (serial), and Linux-mpi
#
#-----------------------------------------------------------------------
#
null  :=
space := $(null) $(null)
comma := $(null),$(null)
#
# Assign defaults:
#
# Determine platform 
ifeq ($(UNAMES),$(null))
  UNAMES := $(shell uname -s)
endif
ifeq ($(EXECNAME),$(null))
  EXECNAME := tgcm
endif
ifeq ($(EXECDIR),$(null))
  EXECDIR := .
endif 
ifeq ($(LIB_NETCDF),$(null))
  LIB_NETCDF := /usr/local/lib
endif
ifeq ($(INC_NETCDF),$(null))
  INC_NETCDF := /usr/local/include
endif
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
# Primary target:
#
all: $(EXECDIR)/$(EXECNAME)

CURDIR   := $(shell pwd)
SOURCES  := $(shell cat Srcfiles)
OBJS     := $(addsuffix .o, $(basename $(SOURCES)))
LIBS     := -L$(LIB_NETCDF) -lnetcdf
LDFLAGS  :=

Depends: Filepath Srcfiles
	perl MkDepends Filepath

Srcfiles: Filepath
	perl MkSrcfiles
#
# Note FFLAGS and CPPFLAGS are included by gmake.
#
$(EXECDIR)/$(EXECNAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

#-----------------------------------------------------------------------
ifeq ($(UNAMES),AIX)

CPPFLAGS   := -WF,-DAIX,-DMSS,-DMPI
OPTIM      := -O3 -qstrict
ifeq ($(DEBUG),TRUE)
  CPPFLAGS := -WF,-DAIX,-DMSS,-DMPI,-DDEBUG
  OPTIM    := -g
  DBGFLAGS := -qinitauto=7FF7FFFF -qflttrap=zero:inv:en -qsigtrap=xl__trcedump
endif
FFLAGS     := $(cpp_path) $(CPPFLAGS) $(OPTIM) $(DBGFLAGS) -qarch=auto -qrealsize=8
FREEFLAGS  := -qsuffix=f=f90:cpp=F90
FIXEDFLAGS := -qfixed
ifeq ($(DEBUG),FALSE)
  LIBS       += -L/usr/local/lib64/r4i4 -lmass
endif
LDFLAGS    += -bloadmap:loadmap -bmaxdata:0X70000000
FC         := mpxlf_r

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

endif
#-----------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)

CPPFLAGS   := -DIRIX -DMSS -DMPI
OPTIM      := -O3
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG
  OPTIM    := -g
endif
FFLAGS     := $(cpp_path) $(OPTIM) -r8 -macro_expand -OPT:Olimit=4479 \
              $(CPPFLAGS) -v
LIBS       += -L/usr/lib32/mips4 -lffio -lcomplib.sgimath -lblas_mp -lfpe
LIBS       += -L/lusr/local/lib32/r4i4 -lmpi
LDFLAGS    += -mips4 -r10000
FC         := f90

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -freeform $<

endif
#-----------------------------------------------------------------------
ifeq ($(UNAMES),Linux)

# 10/04: works at hao under Redhat 9.

CPPFLAGS := -DLINUX -DMSS
OPTIM    := -O2
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG
  OPTIM    := -g
endif
#
# Note: Apparently gmake under Linux includes both FFLAGS and CPPFLAGS by 
#       default, so CPPFLAGS is not included in FFLAGS here.
#
FFLAGS   := $(cpp_path) $(OPTIM) -r8
FC       := pgf90

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -Mfreeform $<

endif
#-----------------------------------------------------------------------
ifeq ($(UNAMES),Linux-mpi)
#
# For SCD linux cluster "lightning"
#
CPPFLAGS := -DLINUX -DMPI -DMSS
# 1/05 btf: fastsse may cause problems in dynamo. Use -fast instead.
#OPTIM    := -O2 -fastsse 
 OPTIM    := -O2 -fast
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG
  OPTIM    := -g
endif
# 5/31/05: path to mpif.h was changed from /opt to /usr after lightning OS upgrade:
#cpp_path += -I/opt/mpich-gm/include
cpp_path += -I/usr/local/mpich-gm/include
#
# Note: Apparently gmake under Linux includes both FFLAGS and CPPFLAGS by 
#       default, so CPPFLAGS is not included in FFLAGS here.
#
FFLAGS   := $(cpp_path) $(OPTIM) -r8
#FC      := pgf90
FC       := mpif90

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -Mfreeform $<

endif
#-----------------------------------------------------------------------
#
RM := rm
clean:
	$(RM) -f Depends Srcfiles *.o *.mod $(EXECDIR)/$(EXECNAME)
#
# Include object file targets and dependencies made by MkDepends:
#
include Depends
