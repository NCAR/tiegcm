#
# OS    = Operating System from uname
# On Input:
#   File Srcfiles contains list of source files
#   File $(MAKE_MACHINE) (e.g. generic_hao64) contains compiler specific Make vars
#     INC_NETCDF: Include directory for netcdf
#     LIB_NETCDF: Library directory for netcdf
#     DBGFLAGS: Debug flags
#   File Make.env contains additional Make vars for the current run:
#     DIRS:     List of source-code directories (in Make.env)
#     DEBUG:    TRUE or FALSE
#
SHELL = /bin/sh

include Make.env         # made by job script before gmake
include $(MAKE_MACHINE)  # machine-specific file with compiler/make info

null     :=
comma    := $(null),$(null)
space    := $(null) $(null)

OS        := $(shell uname -s)
VPATH     := $(subst $(space),:,$(DIRS))
INC_DIRS  := $(foreach dir,$(DIRS),-I$(dir)) -I$(INC_NETCDF)
SOURCES   := $(shell cat Srcfiles)
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))
# LIBS may be specified in Make.machine file, so append here
LIBS      += -L$(LIB_NETCDF) -lnetcdf 
COUPLING  := NO
FC        := $(F90)

all: Depends

Depends::
	@echo "OS           =$(OS)" 
	@echo "MAKE_MACHINE =$(MAKE_MACHINE)"
	@echo "DIRS         =$(DIRS)"
	@echo "INC_DIRS     =$(INC_DIRS)"
	@echo "VPATH        =$(VPATH)"
	@echo "FC           =$(FC)"
	@echo "LIB_NETCDF   =$(LIB_NETCDF)"
	@echo "INC_NETCDF   =$(INC_NETCDF)"
	@echo "EXECNAME     =$(EXECNAME)"
	@echo "DEBUG        =$(DEBUG)"
	perl Mkdepends $(DIRS)
	$(MAKE) $(EXECNAME)

$(EXECNAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

ifeq ($(MPI),TRUE)
  CPPFLAGS += -DMPI
  FC := $(MPIF90)
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
ifeq ($(DEBUG),TRUE)
  CPPFLAGS += -DDEBUG
  OPTIM    := -g
endif

ifeq ($(OS),Linux)
  CPPFLAGS := -DLINUX $(CPPFLAGS)
endif

ifeq ($(OS),AIX)
  CPPFLAGS := -WF -DAIX $(CPPFLAGS)
  CPPFLAGS := $(subst $(space),$(comma),$(CPPFLAGS))
endif

FFLAGS += $(INC_DIRS) $(OPTIM) $(DBGFLAGS)

clean:
	rm -f Depends Srcfiles *.o *.mod $(EXECNAME)
veryclean:
	rm -f Depends Srcfiles *.o *.mod $(EXECNAME) Makefile Make.env $(MAKE_MACHINE) machines.ini Mkdepends defs.h

include Depends

