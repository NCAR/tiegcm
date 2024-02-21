#!/bin/bash
#
# Job script to build and execute TIEGCM on the NCAR supercomputer derecho,
# or similar system running PBS (Load Sharing Facility)
#
# This script should be executed on the derecho command line.
# It will compile in execdir on the interactive login node, then create and
# submit an PBS script based on the #PBS resources specified below.
#
# User should set shell variables and #PBS settings below:
#
#   modeldir:  Root directory to model source (may be an SVN working dir)
#   execdir:   Directory in which to build and execute (will be created if necessary)
#   tgcmdata:  Directory in which startup history and data files are accessed.
#              (If tgcmdata is not set, the model will use env var TGCMDATA)
#   utildir:   Directory containing supporting scripts (default $modeldir/scripts)
#   input:     Namelist input file for the chosen model resolution
#   output:    Stdout file from model execution (will be created)
#   horires:   Model resolution (5, 2.5, 1.25 or 0.625 degrees)
#   vertres:   Model resolution (0.5, 0.25, 0.125 or 0.0625 degrees)
#   zitop:     Pressure level at model upper boundary
#   nres_grid: Number of levels for multigrid solver (magnetic grid)
#   coupling:  Whether to couple with GAMERA via MPI
#   hidra:     Whether to couple with HIDRA via MPI
#   make:      Build file with platform-specific compile parameters (in scripts dir)
#   debug:     If TRUE, build and execute a "debug" run (debug compiler flags are set)
#   execute:   If TRUE, execute the model (build only if execute is FALSE)
#   runscript: PBS script with run commands (submitted with qsub from execdir)
#
# To switch to 2.5x0.25-deg resolution, set horires below to 2.5 and vertres to 0.25,
# and change execdir, tgcmdata and namelist input if necessary.
# Also reset number of processors accordingly below (#PBS -l).
#
modeldir="/glade/campaign/hao/itmodel/tiegcm3.0/tiegcm"
execdir="/glade/derecho/scratch/$USER/tiegcm"
tgcmdata="/glade/campaign/hao/itmodel/tiegcm3.0/data"
utildir="$modeldir/scripts"
input="$utildir/tiegcm_default.inp"
output="tiegcm.out"
horires="2.5"
vertres="0.25"
zitop="7"
mres="2"
coupling="FALSE"
hidra="FALSE"
make="Make.intel_de"
debug="FALSE"
execute="FALSE"
runscript="run.pbs"

project_code="P28100036"
queue="main"
job_priority="regular"
nnodes=1
ncores_per_node=128
walltime="12:00:00"
nprocs=$(($nnodes*$ncores_per_node))

execdir=`perl $utildir/abspath $execdir`
runscript=`perl $utildir/abspath $runscript`
if [ ! -d $execdir ]
then
  mkdir -p $execdir
fi
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                               End user settings
#                        Shell Script for TIEGCM Linux job
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Necessary modules on derecho
module load mkl
module swap netcdf netcdf-mpi
module load esmf
#
mycwd=`pwd`
echo "" ; echo "${0}:"
echo "  Begin execution at `date`"
echo "  Current working directory: $mycwd"
echo "  System: `uname -a`"
echo ""
#
# Verify directories and make_machine file (make execdir if necessary).
# Get absolute path for dirs that are accessed from the execdir.
#
if [ ! -d $modeldir ]
then
  echo ">>> Cannot find model directory $modeldir <<<"
  exit 1
fi
#
# Executable file name:
modelexe=tiegcm.exe

if [ ! -d $utildir ]
then
  echo ">>> Cannot find model scripts directory $utildir <<<"
  exit 1
fi
utildir=`perl $utildir/abspath $utildir`

srcdir=$modeldir/src
if [ ! -d $srcdir ]
then
  echo ">>> Cannot find model source directory $srcdir <<<"
  exit 1
fi
srcdir=`perl $utildir/abspath $srcdir`
#
# Set data directory:
#   If tgcmdata is set above, use it to set (override) env var TGCMDATA
#   If tgcmdata not set above, use env var TGCMDATA
#   If neither are set, print a warning and set both to cwd.
#   Finally, print warning if the directory does not exist.
#
if [ -v tgcmdata ]
then
  tgcmdata=`perl $utildir/abspath $tgcmdata`
  export TGCMDATA=$tgcmdata
  echo "Set env var TGCMDATA = $TGCMDATA"
else
  if [ -v TGCMDATA ]
  then
    tgcmdata=$TGCMDATA
    echo "Set tgcmdata = $TGCMDATA"
  else
    echo "WARNING: TGCMDATA is not set - using cwd"
    tgcmdata=.
    export TGCMDATA=$tgcmdata
  fi
fi
if [ ! -d $tgcmdata ]
then
  echo ">>> Cannot find data directory $tgcmdata"
fi

if [ $horires != "5" ] && [ $horires != "2.5" ] && [ $horires != "1.25" ] && [ $horires != "0.625" ]
then
  echo ">>> Unknown model horizontal resolution $horires <<<"
  exit 1
fi
if [ $vertres != "0.5" ] && [ $vertres != "0.25" ] && [ $vertres != "0.125" ] && [ $vertres != "0.0625" ]
then
  echo ">>> Unknown model vertical resolution $vertres <<<"
  exit 1
fi

if [ $mres == "2" ]
then
  nres_grid=5
elif [ $mres == "1" ]
then
  nres_grid=6
elif [ $mres == "0.5" ]
then
  nres_grid=7
else
  echo ">>> Unsupported magnetic resolution $mres <<<"
  exit 1
fi
#
# Copy make files to execdir if necessary:
#
if [ ! -f $execdir/$make ]
then
  cp $utildir/$make $execdir
fi
if [ ! -f $execdir/Makefile ]
then
  cp $utildir/Makefile $execdir
fi
if [ ! -f $execdir/mkdepends ]
then
  cp $utildir/mkdepends $execdir
fi
#
# Make default namelist input file if not provided by user:
#
if [ ! -f $input ]
then
  echo ">>> Cannot find namelist input file $input <<<"
  exit 1
fi

model=$execdir/$modelexe
input=`perl $utildir/abspath $input`
output=`perl $utildir/abspath $output`
util=`perl $utildir/abspath $utildir`
mklogs=$util/mklogs         # Nov, 2015: mklogs rewritten in python
rmbinchars=$util/rmbinchars # Nov, 2015: remove non-ascii chars from stdout files
#
# Report to stdout:
#
if [ -x `command -v svnversion` ]
then
  svn_version=`svnversion $modeldir`
else
  svn_version="[none]"
fi
#svn_revision="tiegcm3.0" # for svn tag instead of revision number

echo -n "  Model directory:   $modeldir"
echo " (SVN revision $svn_revision)"
echo "  Exec directory:    $execdir"
echo "  Source directory:  $srcdir"
echo "  Data directory:    $tgcmdata"
echo "  Make machine file: $make"
echo "  Namelist input:    $input"
echo "  Stdout Output:     $output"
echo "  Model resolution:  $horires x $vertres"
echo "  Coupling:          $coupling"
echo "  Hidra:             $hidra"
echo "  Debug:             $debug"
#
# If coupling flag has changed from last gmake, clean execdir
# and reset coupling file:
#
if [ -f $execdir/coupling ]
then
  lastcoupling=`cat $execdir/coupling`
  if [ $lastcoupling != $coupling ]
  then
    echo "Clean execdir $execdir because coupling flag switched from $lastcoupling to $coupling"
    mycwd=`pwd` ; cd $execdir ; gmake clean ; cd $mycwd
    echo $coupling > $execdir/coupling
  fi
else
  echo $coupling > $execdir/coupling
  echo "Created file coupling with coupling flag = $coupling"
fi
#
# If hidra flag has changed from last gmake, clean execdir
# and reset hidra file:
#
if [ -f $execdir/hidra ]
then
  lasthidra=`cat $execdir/hidra`
  if [ $lasthidra != $hidra ]
  then
    echo "Clean execdir $execdir because hidra flag switched from $lasthidra to $hidra"
    mycwd=`pwd` ; cd $execdir ; gmake clean ; cd $mycwd
    echo $hidra > $execdir/hidra
  fi
else
  echo $hidra > $execdir/hidra
  echo "Created file hidra with hidra flag = $hidra"
fi
#
# If debug flag has changed from last gmake, clean execdir
# and reset debug file:
#
if [ -f $execdir/debug ]
then
  lastdebug=`cat $execdir/debug`
  if [ $lastdebug != $debug ]
  then
    echo "Clean execdir $execdir because debug flag switched from $lastdebug to $debug"
    mycwd=`pwd` ; cd $execdir ; gmake clean ; cd $mycwd
    echo $debug > $execdir/debug
  fi
else
  echo $debug > $execdir/debug
  echo "Created file debug with debug flag = $debug"
fi
#
# Copy defs header file to execdir, if necessary, according to
# requested resolution. This should seamlessly switch between
# resolutions according to $horires and $vertres.
#
defs=defs.h

cat << EOF > $defs
#define DLAT $horires
#define DLON $horires
#define GLON1 -180
#define DLEV $vertres
#define ZIBOT -7
#define ZITOP $zitop
#define NRES_GRID $nres_grid
EOF

if [ -f $execdir/defs.h ]
then
  cmp -s $execdir/defs.h $defs
  status=$?
  if [ $status -ne 0 ]
  then # files differ -> switch resolutions
    echo "Switching defs.h for model resolution $horires x $vertres"
    mycwd=`pwd` ; cd $execdir ; gmake clean ; cd $mycwd
    cp $defs $execdir/defs.h
  else
    echo "defs.h already set for model resolution $horires x $vertres"
  fi
else # defs.h does not exist in execdir -> copy appropriate defs file
  echo "Copying $defs to $execdir/defs.h for resolution $horires x $vertres"
  cp $defs $execdir/defs.h
fi
#
# cd to execdir and run make:
#
cd $execdir
status=$?
if [ $status -ne 0 ]
then
  echo ">>> Cannot cd to execdir $execdir"
  exit 1
fi
echo ""
echo "Begin building $model in `pwd`"
#
# Build Make.env file in exec dir, containing needed env vars for Makefile:
#
cat << EOF > Make.env
MAKE_MACHINE  = $make
DIRS          = . $srcdir
EXECNAME      = $model
NAMELIST      = $input
OUTPUT        = $output
COUPLING      = $coupling
HIDRA         = $hidra
DEBUG         = $debug
SVN_REVISION  = $svn_revision
EOF
#
# Build the model:
gmake -j8 all
status=$?
if [ $status -ne 0 ]
then
  echo ">>> Error return from gmake all"
  exit 1
fi
#
# Set PBS resource usage (create the runscript in execdir):
# (run commands are appended to this script below)
#
# Set data directory in PBS script:
#
# MPI/PBS job: append mpirun.command command to PBS script
# (it has #PBSs from above)
#
cat << EOF > $runscript
#!/bin/bash
#PBS -N tiegcm
#PBS -A $project_code
#PBS -q $queue
#PBS -j oe
#PBS -l job_priority=$job_priority
#PBS -l select=${nnodes}:ncpus=${ncores_per_node}:mpiprocs=${ncores_per_node}
#PBS -l walltime=$walltime

export TGCMDATA=$tgcmdata

# Execute:
cd $execdir
`cat mpirun.command` -np $nprocs $model $input >& $output

# Save stdout:
$rmbinchars $output # remove any non-ascii chars in stdout file
$mklogs $output     # break stdout into per-task log files
cd `dirname $output`

# Make tar file of task log files:
tar -cf $output.tar *task*.out
echo "Saved stdout tar file $output.tar"
rm *task*.out
EOF
#
if [ $execute == "TRUE" ]
then
  echo " "
  echo "Submitting PBS script $runscript for MPI run"
  qsub < $runscript
else
  echo "Am NOT executing the model because execute = $execute"
fi
#
exit 0
