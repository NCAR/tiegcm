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
# Determine platform 
#
ifeq ($(UNAMES),$(null))
  UNAMES := $(shell uname -s)
endif
#ifeq ($(UNAMEN),$(null))
#  UNAMEN := $(shell uname -n| awk '/^ln[[:digit:]]+en/{print "lightning"}')
#endif
#
# Assign defaults:
#
ifeq ($(EXECNAME),$(null))
  EXECNAME := tgcm.exe
endif
ifeq ($(EXECDIR),$(null))
  EXECDIR := .
endif 
ifeq ($(UTILDIR),$(null))
  UTILDIR := .
endif
#
# Default netcdf lib and include file locations are
# platform dependent:
#
# Netcdf locations on IBM AIX:
ifeq ($(UNAMES),AIX)
 ifeq ($(LIB_NETCDF),$(null))
   LIB_NETCDF := /usr/local/lib32/r4i4
 endif
 ifeq ($(INC_NETCDF),$(null))
   INC_NETCDF := /usr/local/include
 endif
 ifeq ($(MPI),$(null))
   MPI := TRUE
 endif
endif
#
# Netcdf locations on SGI IRIX64:
ifeq ($(UNAMES),IRIX64)
 ifeq ($(LIB_NETCDF),$(null))
   LIB_NETCDF := /usr/local/lib32/r4i4
 endif
 ifeq ($(INC_NETCDF),$(null))
   INC_NETCDF := /usr/local/include
 endif
endif
#
ifeq ($(UNAMES),Linux-mpi)
#
# Netcdf locations on Linux lightning machine:
  ifeq ($(LIB_NETCDF),$(null))   # lightning
    LIB_NETCDF := /usr/local/netcdf/3.5.1-pgi5.2-64/lib
  endif
  ifeq ($(INC_NETCDF),$(null))
    INC_NETCDF := /usr/local/netcdf/3.5.1-pgi5.2-64/include
  endif
  ifeq ($(MPI),$(null))
    MPI := TRUE
  endif
endif
ifeq ($(UNAMES),Linux)
#
# Netcdf locations on hao Linux systems:
  ifeq ($(LIB_NETCDF),$(null))
    LIB_NETCDF := /opt/local/lib
  endif
  ifeq ($(INC_NETCDF),$(null))
    INC_NETCDF := /opt/local/include
  endif
  ifeq ($(MPI),$(null))
    MPI := FALSE
  endif
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
# Default goal (force Mkdepends before making exec):
#
all: Depends 
#
# 9/6/05: MAKEFILE_LIST is empty if bluesky /usr/bin/gmake is used, 
#         but is properly defined by /usr/local/bin/gmake. For now, 
#         I have put full /usr/local/bin path in *ibm.job job scripts.
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
# Note FFLAGS and CPPFLAGS are included by gmake.
#
$(EXECDIR)/$(EXECNAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

#-----------------------------------------------------------------------
ifeq ($(UNAMES),AIX)

ifeq ($(MPI),TRUE)
  CPPFLAGS   := -WF,-DAIX,-DMSS,-DMPI
else
  CPPFLAGS   := -WF,-DAIX,-DMSS
endif
OPTIM      := -O3 -qstrict
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG,-DEXPO
  CPPFLAGS := $(subst $(space),$(comma),$(CPPFLAGS))
  OPTIM    := -g
  DBGFLAGS := -qinitauto=7FF7FFFF -qflttrap=zero:inv:en -qsigtrap=xl__trcedump
# DBGFLAGS := -C -qzerosize -qinitauto=7FF7FFFF -qflttrap=zero:inv:en \
#             -qsigtrap=xl__trcedump
endif
#
# 6/23/06 btf: Now using -q64 (64-bit addressing), s.a. LDFLAGS below.
#              (double res timegcm needs 64-bit)
#
FFLAGS     := $(cpp_path) $(CPPFLAGS) $(OPTIM) $(DBGFLAGS) \
              -qarch=auto -qrealsize=8 -q64
#
FREEFLAGS  := -qsuffix=f=f90:cpp=F90
FIXEDFLAGS := -qfixed
#
ifeq ($(DEBUG),FALSE)
  LIBS     += -lmass
else
  LIBS     +=          # no -lmass if debug
endif
#LDFLAGS    += -bloadmap:loadmap -bmaxdata:0X70000000
 LDFLAGS    += -bloadmap:loadmap -q64
ifeq ($(MPI),TRUE)
  FC         := mpxlf_r
else
  FC         := xlf90
endif

.SUFFIXES:
.SUFFIXES: .F .F90 .c .o

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<
#
# -qsave is necessary for pdedif (glbmean model only):
#    Sets the default storage class for	local variables	to STATIC
#
pdedif.o: pdedif.F
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) -qsave $<

endif
#-----------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)

ifeq ($(MPI),TRUE)
  CPPFLAGS   := -DIRIX -DMSS -DMPI
else
  CPPFLAGS   := -DIRIX -DMSS
endif
OPTIM      := -O3
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG,-DEXPO
  OPTIM    := -g
endif
FFLAGS     := $(cpp_path) $(OPTIM) -r8 -macro_expand -OPT:Olimit=4479 \
              $(CPPFLAGS) -v
LIBS       += -L/usr/lib32/mips4 -lffio -lcomplib.sgimath -lblas_mp -lfpe
ifeq ($(MPI),TRUE)
  LIBS       += -L/lusr/local/lib32/r4i4 -lmpi
endif
#
LDFLAGS    += -mips4 -r10000
FC         := f90

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -freeform $<
#
# -static is necessary for pdedif (glbmean model only):
#    Sets the default storage class for	local variables	to STATIC
#
pdedif.o: pdedif.F
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) -static $<

endif
#-----------------------------------------------------------------------
ifeq ($(UNAMES),Linux-mpi)
#
# For SCD linux cluster "lightning"
#
ifeq ($(MPI),TRUE)
  CPPFLAGS := -DLINUX -DMPI -DMSS
#
# 5/31/05: path to mpif.h was changed from /opt to /usr after lightning OS upgrade:
# cpp_path += -I/opt/mpich-gm/include
  cpp_path += -I/usr/local/mpich-gm/include
else
  CPPFLAGS := -DLINUX -DMSS
endif

# 1/05 btf: fastsse may cause problems in dynamo. Use -fast instead.
#OPTIM    := -O2 -fastsse 
 OPTIM    := -O2 -fast
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG,-DEXPO
  OPTIM    := -g
endif
#
# Note: Apparently gmake under Linux includes both FFLAGS and CPPFLAGS by 
#       default, so CPPFLAGS is not included in FFLAGS here.
#
FFLAGS   := $(cpp_path) $(OPTIM) -r8
#
 LIBS    += 
#LDFLAGS +=

ifeq ($(MPI),TRUE)
  FC := mpif90
else
  FC := pgf90
endif

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -Mfreeform $<
#
# -Msave is necessary for pdedif (glbmean model only):
#    Sets the default storage class for	local variables	to STATIC
#
pdedif.o: pdedif.F
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) -Msave $<
endif
#
#-----------------------------------------------------------------------
ifeq ($(UNAMES),Linux)
#
# Non-lightning Linux (probably hao Linux):
#
# 6/05: This works for serial linux with Fedora core 3 (sundog, suncat, etc)
#       Trying to make it work w/ MPI, e.g., on hao 4-proc machine "tor".

CPPFLAGS := -DLINUX -DMSS
ifeq ($(MPI),TRUE)
  CPPFLAGS := -DLINUX -DMPI -DMSS
  cpp_path += -I/usr/local/include
endif
OPTIM    := -O2
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG,-DEXPO
  OPTIM    := -g
endif
#
# Note: Apparently gmake under Linux includes both FFLAGS and CPPFLAGS by 
#       default, so CPPFLAGS is not included in FFLAGS here.
#
FFLAGS   := $(cpp_path) $(OPTIM) -r8
#
ifeq ($(MPI),TRUE)
  LIBS   += -L/usr/local/lib -lmpich
endif
#LDFLAGS +=
FC       := pgf90

gw_drag.o: gw_drag.F
	$(FC) -c $(FFLAGS) -Mfreeform $<
#
# -Msave is necessary for pdedif (glbmean model only):
#    Sets the default storage class for	local variables	to STATIC
#
pdedif.o: pdedif.F
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) -Msave $<

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

