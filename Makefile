#
UNAMES = $(shell uname -s)
UNAMEM = $(findstring CRAY,$(shell uname -m))
EXEC = tiegcm1
#
#----------------------------------------------------------------------
# CRAY under UNICOS (J90 or J90se):
# 
ifeq ($(UNAMEM),CRAY)
#ifeq ($(UNAME),unicos)
#CRAY_begin
#
FC      = f90
#
FFLAGS = -xmic,dir -m4 -O2 -F overindex,nofastint -I/usr/local/include \
         -DUNICOS -DMSS
LIBS = -L /usr/local/lib -l ncarg,ncarg_gks,ncarg_c,X11,ncarm, \
       -l ncaro,mss,alfpack,ecmfft,fitpack,fishpak,fftpack,hpf,netcdf
LDFLAGS = -Wl"-f zeros -M,stat"
#CRAY_end
endif
#
#----------------------------------------------------------------------
# SGI O2K running IRIX64:
#
# For ncar chinook SGI Origin 3800: 
#
ifeq ($(UNAMES),IRIX64)
FC      = f90
FFLAGS = -mips4 -r8 -O3 -macro_expand -I/usr/local/include \
         -DMSS -DIRIX -DMPI
LIBS = \
  -L /usr/lib32/mips4 -lffio -lcomplib.sgimath -lblas_mp -lfpe \
  -L /usr/local/lib32/r4i4 -lmss -lnetcdf -lmpi \
  -L /usr/local/dcs-3.3/lib -ldcs
LDFLAGS = -mips4 -r10000 -O3
endif
#----------------------------------------------------------------------
# IBM SP: AIX
#
ifeq ($(UNAMES),AIX)
#
# For profiling:
# Use -pg for gprof or xprofiler (must also be in LDFLAGS below)
#   Add -g to -pg for profiling at source level.
#   After execution, gmon.out.x files will be produced. To run
#   gprof on task 0, use: "gprof tgcm.aix gmon.out.0"
#
# Use mpxlf_r for MPI, xlf90_r for OpenMP.
#
# 5/00: use -lmass with care -- it can cause fpe trapping problems.
# 12/4/01: -lmass works with tiegcm1 if -qflttrap and -qsigtrap are
#          removed. The MASS lib greatly improves performance.
# 3/02: Since upgrading to xlf7.1.1.2 consultants recommend avoiding
#       -qtune and -qarch (which were both set to pwr3).
#
# - - - - - - - - - - - - - - - - - 
#
# For MPI:
#
 FC = mpxlf_r
#
# "Production mode": Model will not crash with these settings, 
#   but it may produce NaNQ or other invalids which may end 
#   up on the histories. These settings produce the fastest 
#   but least "protected" code.
#
 FFLAGS= -O3 -qstrict -qfixed -qrealsize=8 \
 	 -WF,-DMPI,-DAIX,-DMSS -I/usr/local/include
#
# With these settings, all automatic variables are init to NaNS. 
# The program will stop and core dump with "invalid" fpe on referencing 
#   uninitialized variables or other invalids (e.g., exp()).
# (add -g and remove -qstrict for source code reference debugging)
#
#FFLAGS= -O3 -qstrict -qinitauto=7FF7FFFF -qfloat=nans \
#        -qflttrap=invalid:enable  -qsigtrap=xl__trcedump \
#        -qfixed -qrealsize=8 -WF,-DMPI,-DAIX,-DMSS \
#        -I/usr/local/include
#
# Trapping and debugging: this is 4-5x slower than with -O3 above. 
# (cat set check_exp in cons.F to get around invalid exp() fpe's)
#
#FFLAGS= -g -qinitauto=7FF7FFFF -qfloat=nans \
#        -qflttrap=invalid:enable  -qsigtrap=xl__trcedump \
#        -qfixed -qrealsize=8 -WF,-DMPI,-DAIX,-DMSS \
#        -I/usr/local/include
#
#FFLAGS= -g -pg -O3 -qstrict -qfixed -qrealsize=8 \
#	 -WF,-DMPI,-DAIX,-DMSS -I/usr/local/include
#
 LIBS = -lmass -L /usr/local/lib32/r4i4 -lnetcdf -lmss \
        -L /usr/local/dcs-3.3/lib -ldcs
# - - - - - - - - - - - - - - - - - 
# For MPI with VAMPIRtrace:
# Put the following in .cshrc:
# setenv PAL_ROOT /usr/local/vampir
# setenv PAL_LICENSEFILE /usr/local/vampir/etc/license.dat
#
#FC = mpxlf_r
#VTroot = /usr/local/vampir
#FFLAGS= -qfixed -qrealsize=8 \
#        -O3 -qstrict -WF,-DMPI,-DAIX,-DMSS,-DVT \
#        -I/usr/local/include -I$(VTroot)/include
#LIBS = -lmass -L /usr/local/lib32/r4i4 -lnetcdf -lmss \
#       -L /usr/local/dcs-3.2/lib -ldcs \
#       -L$(VTroot)/lib -lVT -lld
# - - - - - - - - - - - - - - - - - 
# For MPI and OpenMP:
#FC = mpxlf_r
#FFLAGS= -qfixed -qrealsize=8 \
#        -qsmp=omp -qthreaded \
#	 -O3 -qstrict -WF,-DMPI,-DAIX,-DMSS \
#	 -qflttrap=zero:invalid:enable -qsigtrap=xl__trcedump \
# -I/usr/local/include
#LIBS = -L /usr/local/lib32/r4i4 -lnetcdf -lmss -lxlsmp
#LDFLAGS = -qsmp=omp
# - - - - - - - - - - - - - - - - - 
#
# For OpenMP: 
# 4/19/00: -O5[4] -qstrict -qsmp=omp -qthreaded compiles, but 
#   crashes at end of lamdas. Reducing to -O3 -qstrict works.
# 6/9/00: is not crashing, but produces garbage NaN's
#   (is it necessary to make sure 4 threads are all on same node?)
#
#FC = xlf90_r
#FFLAGS= -qfixed -qrealsize=8 \
#	-qsmp=omp -qthreaded \
#	-qalias=noaryovrlp -O3 -qstrict -WF,-DAIX,-DMSS \
#	-qflttrap=zero:enable -qsigtrap=xl__trcedump \
#	-I/usr/local/include
#LIBS    = -L /usr/local/lib32/r4i4 -lnetcdf -lmss -lxlsmp
# - - - - - - - - - - - - - - - - - 
#
# For serial (no OpenMP or MPI) debug:
#
#FC = xlf90
#FFLAGS= -qfixed -qrealsize=8 \
#	-qalias=noaryovrlp -g -WF,-DAIX,-DMSS \
#	-qflttrap=zero:invalid:enable -qsigtrap=xl__trcedump \
#	-I/usr/local/include
#LIBS = -L /usr/local/lib32/r4i4 -lnetcdf -lmss
# - - - - - - - - - - - - - - - - - 
#
# -bmaxdata:0X70000000 is necessary for dz=0.25 vertical resolution.
 LDFLAGS = -bloadmap:loadmap -bmaxdata:0X70000000
#LDFLAGS = -pg -bloadmap:loadmap -bmaxdata:0X70000000
endif
#----------------------------------------------------------------------
# COMPAQ: OSF1 prospect (8 nodes, 4 pe/node)
# Use gmake. As of 6/00, no batch system is available (interactive only)
#
ifeq ($(UNAMES),OSF1)
#OSF1_begin
FC      = f90
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For threaded OpenMP:
#
# 3/15/00: prospect generates internal compiler error if -omp,
#          -O5, and -fast are set. Works if -fast is removed.
# 4/00: works only with -O3 (and no -fast)
#
#FFLAGS  = -arch ev6 -omp -real_size 64 -O3 \
#	-I/usr/local/include -DOSF1 -DMSS -align dcommons
#LDFLAGS = -omp
#LIBS    = -L/usr/local/lib -lnetcdf -lmss -lelan
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For mpi only:
# e.g., for 12 mpi tasks (no omp): ncpus=12, ntasks (nodes) = 3
# prun -n $ncpus -N $ntasks -t -v sh -c "$exec < $inp"
#
# NCAR's prospect:
#FFLAGS  = -arch ev6 -real_size 64 -O3 \
#  -I/usr/local/include -DOSF1 -DMPI -DMSS -align dcommons
#
# PSC's lemieux (note no -DMSS):
# Must execute "module load netcdf" so set the NETCDF env vars
#   (this can be in .cshrc and/or the batch job script)
 FFLAGS  = -tune host -arch host -real_size 64 -O3 -DOSF1 -DMPI -align dcommons \
   -I$(NETCDF_INC) -I/usr/include
 LIBS    = -lmpi -lelan -L$(NETCDF_LIB) -lnetcdf
#LIBS    = -L/usr/local/lib -lmpi -lelan -L$(NETCDF_LIB) -lnetcdf
#LDFLAGS = -lmpi -lelan -lnetcdf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For hybrid MPI/OpenMP:
#
# Batch PBS jobs using prun command:
# OMP_NUM_THREADS should be set to the same number as -c in prun command.
# For example, to make a 12-proc (3 nodes 4 pe's/node) run:
# prun -n3 -c4 -t -v sh -c "tgcm.osf < tgcm14a.inp"
# (Use rinfo -a to see current activity on the machine)
#
# Can use -fpe2 to trap fpe's (e.g. overflows in chapmn and underflows in qrj),
#   but is slow. Can use -fpe0 -synchronous_exceptions to debug fpe's.
#
#FFLAGS  = -arch ev6 -omp -real_size 64 -O3 -fpe0 -synchronous_exceptions \
#	-I/usr/local/include -DOSF1 -DMPI -DMSS -align dcommons
#FFLAGS  = -arch ev6 -omp -real_size 64 -O3 \
#-I/usr/local/include -DOSF1 -DMPI -DMSS -align dcommons
#LDFLAGS = -omp -lmpi -lelan
#LIBS    = -L/usr/local/lib -lnetcdf -lmss -lmpi -lelan
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For single-task runs for interactive debug (no omp or mpi):
# (interactive multi-task prun is not allowed) 
# Do not use fixed object code -- recompile entire code
#   (set COMPILE_ALL in run script)
# The default -fpe0 aborts on fpe's. Use -fpe2 to trap fpe's.
#
#FFLAGS  = -arch ev6 -real_size 64 -fpe0 -O3 \
#	-I/usr/local/include -DOSF1 -DMSS -align dcommons
#LDFLAGS =
#LIBS    = -L/usr/local/lib -lnetcdf -lmss -lelan
#
#OSF1_end
endif
#----------------------------------------------------------------------
# Sun:
#
ifeq ($(UNAMES),SunOS)
#SunOS_begin
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For scd k2 Sun:
# fpversion shell command says -xtarget=ultra2 can be used.
#FC      = /opt/fsc/SUNWspro/bin/f90
#
# For serial:
#FFLAGS  = -g -O -xtypemap=real:64,integer:32 -I/fs/local/include -DSUN
#LDFLAGS = 
#
# For OpenMP:
#FFLAGS  = -O3 -xtypemap=real:64,integer:32 -I/fs/local/include -DSUN \
#          -xtarget=ultra2 -ftrap=no%overflow -DMSS \
#          -mp=openmp -explicitpar -stackvar
#LDFLAGS = -O3 -explicitpar -mp=openmp
#LIBS    = -L/fs/local/lib -lnetcdf
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For hao Sun:
#
 FC      = /opt/SUNWspro/bin/f90
 FFLAGS  = -g -I/opt/share/PUBLIC/include -xpp=cpp -DSUN -DMSS -xtypemap=real:64
#FFLAGS  = -O -I/opt/share/PUBLIC/include -xpp=cpp -DSUN -DMSS -xtypemap=real:64
#LIBS    = -L/opt/share/PUBLIC/lib.SunOS -lnetcdf3 -lm \
#          -L/opt/local/lib -ldcs
#
# This LDLIBS worked on Solaris 5.8 with Huixin's code, which includes 
# apex code:
#IDLDIR=/opt/share/idl
#BINDIR=$(IDLDIR)/bin/bin.solaris2.sparc
#LDFLAGS = -xlic_lib=sunperf -L$(BINDIR) \
#   -lidl -lXm -lXext -lXt -lX11 -lcurses -lF77
#
# Load Sun Performance Lib with "-xlic_lib=sunperf" -- this has 
# Sun optimized versions of lapack, blas, linpack, fftpack.
#
 LIBS    = -L/opt/share/PUBLIC/lib.SunOS -lnetcdf3 -lm \
           -L/opt/local/lib -ldcs -xlic_lib=sunperf
#SunOS_end
endif
#----------------------------------------------------------------------
# Linux:
#
ifeq ($(UNAMES),Linux)
 FC = pgf90
#
# Redhat 8 needs -Msecond_underscore, but may get undefined __ctype_b
#FFLAGS = -I/opt/local/include -Msecond_underscore
#LIBS = -L/opt/local/lib -lnetcdf
#
# Redhat 9 works (do not use -Msecond_underscore):
 FFLAGS = -r8 -O2 -I/opt/local/include -DMSS -DLINUX
#FFLAGS = -g -r8 -I/opt/local/include -DMSS -DLINUX
 LIBS = -L/opt/local/rh9/lib -lnetcdf -L /opt/local/dcs-3.3/rh9/lib -ldcs
endif
#----------------------------------------------------------------------
#EXEC_begin
#
include Objects
#OBJS = $(FSRCS:.f=.o) $(FSRCS:.F=.o)

.SUFFIXES: .F

$(EXEC):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

include Depends

.f.o:
	$(FC) -c $(FFLAGS) $<
.F.o:
	$(FC) -c $(FFLAGS) $<
#
# flint is available on babyblue as of 1/02, path = /home/babyblue/ipt/fl437
# (use "startlm_flint" to start flint license manager)
#
flint_basic:
	flint -I . -I /usr/local/include -9 -Mdepend -u -g *.F
flint_tree:
	flint -I . -I /usr/local/include -9 -Mdepend -Ttrim,condensed *.F
flint_xref:
	flint -I . -I /usr/local/include -9 -Mdepend -x -g *.F
flint_all:
	flint -I . -I /usr/local/include -9 -Mdepend -x -t -u -g *.F
#EXEC_end
