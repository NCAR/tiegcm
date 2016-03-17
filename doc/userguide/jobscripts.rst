
.. _jobscript:

Using the job scripts to set up and submit a model run 
======================================================


The Linux desktop job script (tiegcm-linux.job)
-----------------------------------------------

Take a look at the default Linux job script 
:download:`tiegcm-linux.job <../../scripts/tiegcm-linux.job>`. Near the top are
several shell variables, with their default settings, which configure the 
job script (variables and values may vary between model versions)::

  set modeldir = tiegcm_trunk
  set execdir  = /hao/aim/$user/tiegcm_trunk/tiegcm.exec
  set tgcmdata = /hao/aim/tgcm/data/tiegcm2.0
  set input    = $modeldir/scripts/tiegcm_res5.0_default.inp
  set output   = tiegcm.out
  set make     = Make.intel_hao64
  set modelres = 5.0
  set mpi      = TRUE  # must be TRUE for tiegcm2.0 and later
  set nproc    = 4
  set debug    = FALSE
  set exec     = TRUE
  set utildir  = $modeldir/scripts

Following are brief explanations of the job script shell variables:

.. note::
   
   Absolute or relative paths are acceptable when specifying directories.
   Relative paths should be relative to the *working directory* (:term:`workdir`).
   In practice, modeldir is usually relative to the working directory, and
   execdir and tgcmdata are usually absolute paths.

.. describe:: modeldir

   The model root directory (:term:`modeldir` from the source code download). 
   The example above assumes the user has checked out the trunk revision as "tiegcm_trunk".
   This directory contains subdirectories :term:`src/` , :term:`scripts/` , :term:`doc/` , 
   :term:`tgcmrun/`, and :term:`benchmarks/`.

.. describe:: execdir

   This is the execution directory (:term:`execdir`), in which the model will
   be built and executed. It will be created if it does not already exist.
   It is typically on a large temporary disk. This directory will also contain 
   the model output :term:`netCDF` history files (see also :ref:`historyoutput`)

.. describe:: tgcmdata

   Directory containing startup history files and data files for model input.
   It is normally on a large temporary disk.  These files are available from 
   the :ref:`data download tar file <download>` (separate downloads for each
   model resolution).  Note that setting :term:`tgcmdata` in the job script is 
   optional: if it is specified, it will override any setting of the 
   :term:`TGCMDATA` environment variable.  If it is not specified, the job 
   script will use the TGCMDATA environment variable. If neither are set, 
   tgcmdata will default to the current working directory.

.. describe:: input

   The :ref:`namelist input file <namelist>`. The default namelist file is in the scripts directory
   under the model root with file name :download:`tiegcm_res5.0_default.inp <../../scripts/tiegcm_res5.0_default.inp>`
   (for 5-degree resolution), or :download:`tiegcm_res2.5_default.inp <../../scripts/tiegcm_res2.5_default.inp>`
   (for 2.5-degree resolution).
   The default input file can be copied to the working directory, modified, and renamed for your
   own runs. In that case, be sure to reset the input file in the job script.

.. describe:: make

   Make file containing platform-specific compiler flags, library locations, etc.
   If not otherwise specified with a path, the job script will look for this
   file in the :term:`scripts/` directory. This file is included in the main Makefile
   (scripts/Makefile).  The user can either make necessary adjustments to 
   an existing make file, or write their own for a different platform/compiler 
   system.

   There are three such makefiles available in the :term:`scripts/` directory for the
   Linux desktop platform:

     * :download:`Make.intel_hao64 <../../scripts/Make.intel_hao64>` (for Intel compiler)
     * :download:`Make.pgi_hao64 <../../scripts/Make.pgi_hao64>` (for PGI compiler)
     * :download:`Make.gfort_hao64 <../../scripts/Make.gfort_hao64>` (for gfortran compiler)

   You will need to set the paths to your local netCDF and :term:`ESMF` libraries
   in these makefiles.

.. describe:: output

   Name of the file to receive stdout :term:`output` from the model. If this 
   pre-exists, it will be overwritten when the model is executed.
   Here is an example stdout file from the root mpi task of a 4-processor
   run (5-degree resolution) on a Linux desktop machine:
   :download:`tiegcm_task0000.out <_static/tiegcm_task0000.out>`

.. describe:: mpi

   Logical flag indicating whether or not to link the MPI library for a 
   multi-processor parallel run.

.. warning::

   For |model| versions |version| and later, non-MPI runs (mpi=FALSE) are NOT supported.
   However, mpi runs (mpi=TRUE) with a single processor (nproc=1) ARE supported.

.. describe:: nproc

   Number of processors to use in parallel execution. This will be the number
   of MPI tasks made available for the domain decomposition. On linux desktops,
   this is typically 4. For |model| on linux supercomputer clusters (e.g., the 
   NCAR |ncarsuper| system, where there are 16 processors per node), the recommended 
   number is 16 for 5.0-degree resolution, or 64 for 2.5-degree resolution. 
   For debug purposes, nproc=1 is supported. The models have been tested with
   the following processor counts: 1,4,8,12,16,24,32,48,64,72,80. See 
   :download:`performance table <_static/perf.table>` for performance estimates
   at recommended processor counts and timesteps.

.. _modelres:
.. describe:: modelres

   Model resolution. Two resolutions are supported: 
     * modelres = 5.0 sets 5-degree lat x lon horizontal, and dz=0.50 vertical
     * modelres = 2.5 sets 2.5-degree lat x lon horizontal, and dz=0.25 vertical

   If the resolution is changed, the model should be recompiled before re-executing 
   the job script (type "*gmake clean*" in the :term:`execdir`.

   For more information, see :ref:`Grid Structure and Resolution <resolution>`.

.. describe:: debug
   
   If debug = TRUE, the job script will compile the build with debug flags set.
   Debug flags specific to the compiler are set in the make file. If debug
   is changed, the code should be recompiled (type "gmake clean" in the :term:`execdir`
   before re-executing the job script).

.. describe:: exec

   If exec = TRUE, the job script will execute the model after compilation,
   otherwise, the job script will stop after compilation without execution.

.. describe:: utildir

   The utility directory containing supporting scripts. This is normally the :term:`scripts/`
   subdirectory in the model root directory :term:`modeldir`

You are now ready to build and execute a default run. To do this, simply execute the job script
as follows::

  $ tiegcm-linux.job &


The compilation output will be displayed. If the build is successful (and exec=TRUE),
the model will be executed, and stdout will go to the specified :term:`output` file.
If the job is successful, you can edit and rename the namelist input file, reset
:term:`namelist input file <namelist input>` in the job script, and re-execute the job script. 
If there has been no change to the source code, it will not need to recompile, and will 
use the pre-existing executable.

The |ncarsuper| supercomputer job script (tiegcm-ys.job)
--------------------------------------------------------

.. note::

   This section contains information that is specific to user's
   of the NCAR Linux Supercomputer |ncarsuper|:

   | Linux yslogin3 2.6.32-358.el6.x86_64 #1 SMP 
   | Tue Jan 29 11:47:41 EST 2013 x86_64 x86_64 x86_64 GNU/Linux

   For more information about the NCAR |ncarsuper| system, see
   |ncarsuper_url|

The model can be built and executed on |ncarsuper| using the Intel
compiler and the intelmpi implementation. To do this, copy and modify
the job script tiegcm-ys.job from the scripts directory.

The |ncarsuper| job script :download:`tiegcm-ys.job <../../scripts/tiegcm-ys.job>` 
has the same user-settable shell variables as the Linux job script, but 
the default settings are slightly different::

  set modeldir = tiegcm_trunk
  set execdir  = /glade/scratch/$user/tiegcm_trunk/tiegcm.exec
  set tgcmdata = /glade/p/hao/tgcm/data/tiegcm2.0
  set input    = $modeldir/scripts/tiegcm_res5.0_default.inp
  set output   = tiegcm.out
  set modelres = 5.0
  set make     = Make.intel_ys
  set mpi       = TRUE   # must be TRUE for tiegcm2.0 and later
  set debug     = FALSE    
  set exec      = TRUE
  set utildir   = $modeldir/scripts
  set runscript = run.lsf

In this example, it is assumed the user has checked-out the trunk revision
as "tiegcm_trunk".  Note the :term:`execdir` name, and the ``make`` file 
:download:`Make.intel_ys <../../scripts/Make.intel_ys>`.
The model :term:`resolution` in this case is 5.0 degrees.

Also note the special "#BSUB" directives near the top of the |ncarsuper|
job script (descriptions in the right-hand column are for this document
only, and are not in the script itself)::

  #BSUB -J tiegcm                # job name
  #BSUB -P P28100036             # authorized project number
  #BSUB -q premium               # premium queue
  #BSUB -o tiegcm.%J.out         # stdout file
  #BSUB -e tiegcm.%J.out         # stderr file
  #BSUB -N
  #BSUB -u $LOGNAME@ucar.edu     # email notification address
  #BSUB -W 1:00                  # wallclock limit hours:minutes
  #BSUB -n 16                    # number of processors (mpi tasks)
  #BSUB -R "span[ptile=16]"      # use 16 processors per node

These are resource settings for the Load Sharing Facility (LSF),
the batch queuing system sold by Platform Computing. The LSF is
used for scheduling jobs on the |ncarsuper| system at NCAR.
This job will be submitted to the premium queue command, 
requesting 16 processors with a wallclock limit of 1 hour. 

To submit the |ncarsuper| job, simply execute the job script 
on the command line. It will build the model on the interactive
node, and if successful, the runscript (run.lsf by default) will 
be created and submitted to the LSF via the bsub command. 

Watch the progress of your LSF job with the command::

  $ bjobs

You can kill a LSF job with this command::

  $ bkill job_ID

Where ``job_ID`` is the job identifier given in the ``bjobs`` command.

For more information about the LSF, see the Wikipedia site:

http://en.wikipedia.org/wiki/Platform_LSF

or the Platform Computing site:

http://www.platform.com/workload-management/high-performance-computing/lp

