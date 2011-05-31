
.. _quickstart:

QuickStart Procedure
====================

This document is intended to provide a tutorial-like procedure for downloading 
source code and start-up data files, building the model, and executing a short 
default run on a 64-bit Linux system.

.. note::

   This document is up to date for version |version| of the TIEGCM

.. index:: download

.. _download:

Downloading the model source code and required data files
---------------------------------------------------------

The model source code and related input data files may be downloaded from 
the TIEGCM download page of the main TGCM website:

http://www.hao.ucar.edu/modeling/tgcm/download.php

You will need to provide an email address (login and password are NOT required).
Documentation and Postprocessor codes are also available on the download site,
but all you need for now is the source code, and corresponding data files.
Both are provided as gzipped tar files.

After downloading the two gzipped tar files to an empty working directory (:term:`workdir`) 
on a large disk system (for example, ``/mydisk/tiegcm``), uncompress and extract 
the source code and related scripts and documentation::

    $ gunzip tiegcmx.xx.tar.gz
    $ tar xvf tiegcmx.xx.tar

where ``x.xx`` is the model version downloaded. 

Next, make a directory to hold
the data files (for example, ``/mydisk/tiegcm/data``), and uncompress and extract 
the data tar file into that directory.  Then set environment variable ``$TGCMDATA``,
e.g., for the c-shell, add this line to your .cshrc file::

    setenv TGCMDATA /mydisk/tiegcm/data

At this point, you should have something like this in your working directory
``/mydisk/tiegcm``::

  total 2376
  -rw-rw-r-- 1 user tgcm    4928 Jun  1  2010 README.download
  -rw-r--r-- 1 user tgcm    2886 Jun  1  2010 Release_Notes
  drwxrwxr-x 2 user tgcm    4096 Apr 22 08:52 data/
  -rwxrwxr-x 1 user tgcm   10137 May 31  2010 tiegcm-ibm.job*
  -rwxrwxr-x 1 user tgcm   10671 May 31  2010 tiegcm-linux.job*
  drwxrwxr-x 5 user tgcm    4096 Jun  1  2010 tiegcmx.xx/
  -rw-r--r-- 1 user tgcm    6116 May 31  2010 tiegcmlicense.txt
 
These files and directories contain the following:

.. describe:: README.download

   :download:`README.download <_static/README.download>`
   Instructions for building and making a short default run.

.. describe:: Release_Notes

   Release notes for this version of the model.

.. describe:: data/

   Directory containing the downloaded data files (this is ``$TGCMDATA``)

.. describe:: tiegcm-ibm.job

   Job script for building and executing under IBM/AIX systems.
   (:download:`default tiegcm-ibm.job <_static/tiegcm-ibm.job>`)

   Read more about :ref:`Running the model on IBM/AIX systems <ibm-systems>`.

.. describe:: tiegcm-linux.job 

   Job script for building and executing under Linux (64-bit) systems.
   (:download:`default tiegcm-linux.job <_static/tiegcm-linux.job>`)

.. describe:: tiegcmx.xx/

   Model root directory, containing source code, supporting scripts,
   and documentation.

.. tiegcmlicense.txt 

   Academic license agreement.
   (:download:`tiegcmlicense.txt <_static/tiegcmlicense.txt>`)

You are now prepared to build the model and make a short default run
using the job script.

.. index:: default ; Linux run

.. _jobscript:

Making a default run on a 64-bit Linux system
---------------------------------------------

Take a look at the Linux job script 
:download:`tiegcm-linux.job <_static/tiegcm-linux.job>`. Near the top are
several shell variables, with their default settings, which configure the 
job script (variables and values may vary somewhat between model versions)::

 set modeldir = tiegcmx.xx
 set execdir  = tiegcm-linux
  set make     = Make.intel_hao64
 #set make     = Make.pgi_hao64
 #set input    = tiegcm.inp
 set output   = tiegcm.out
 set mpi      = TRUE
 set nproc    = 4
 set modelres = 5.0
 set debug    = FALSE
 set exec     = TRUE
 set utildir  = $modeldir/scripts

Following are brief explanations of the job script shell variables:

.. note::
   
   Absolute or relative paths are acceptable when specifying directories.
   Relative paths should be relative to the *working directory* (:term:`workdir`).

.. index:: modeldir
.. describe:: modeldir

   The model root directory (:term:`modeldir` from the source code download). 
   This will contain subdirectories ``src/``, ``scripts/``, ``doc/``, etc.

.. index:: execdir
.. describe:: tiegcm-linux

   This is the execution directory (:term:`execdir`), in which the model will
   be built and executed. It will be created if it does not already exist.
   This directory will also contain the model output netCDF history files.

.. index:: make
.. describe:: make

   Make file containing platform-specific compiler flags, library locations, etc.
   If not otherwise specified with a path, the job script will look for this
   file in ``modeldir/scripts``. This file is included in the main Makefile
   (``scripts/Makefile``).  The user can either make necessary adjustments to 
   an existing ``make`` file, or write their own for a different platform/compiler 
   system.

   Here is an example ``make`` file for 64-bit HAO Linux systems using the ifort
   Intel compiler: :download:`Make.intel_hao64 <_static/Make.intel_hao64>`

.. describe:: input

   The namelist input file. When this is commented (as above), the job script
   will make a default namelist file :download:`tiegcm_default.inp <_static/tiegcm_default.inp>`, 
   and use it for the default run. Later, you can edit this file for your own runs, 
   rename it, and reset and uncomment the ``input`` shell variable in the job script.

.. describe:: output

   Name of the file to receive stdout output from the model. If this pre-exists, 
   it will be overwritten when the model is executed.
   
   Here is an example stdout file from a single-processor default run:
   :download:`tiegcm_default.out <_static/tiegcm_default.out>`

.. describe:: mpi

   Logical flag indicating whether or not to link the MPI library for a 
   multi-processor parallel run. If FALSE, the MPI library is not linked,
   and it is assumed the model will be run in serial (single-processor) mode.

.. describe:: nproc

   Number of processors to use in a parallel execution. This is ignored if
   ``mpi`` is FALSE.

.. _modelres:

.. describe:: modelres

   Model resolution. Two resolutions are supported: 
     * modelres = 5.0 sets 5-degree lat x lon horizontal, and dz=0.50 vertical
     * modelres = 2.5 sets 2.5-degree lat x lon horizontal, and dz=0.25 vertical

   If the resolution is changed, the model should be recompiled before re-executing 
   the job script (type "*gmake clean*" in the execdir).

   For more information, see :ref:`Grid Structure and Resolution <resolution>`.

.. describe:: debug
   
   If ``debug`` = TRUE, the job script will compile the build with debug flags set.
   Debug flags specific to the compiler are set in the ``make`` file. If ``debug`` 
   is changed, the code should be recompiled (type "gmake clean" in the ``execdir``
   before re-executing the job script).

.. describe:: exec

   If ``exec`` = TRUE, the job script will execute the model after compilation,
   otherwise, the job script will stop after compilation without execution.

.. describe:: utildir

   The utility directory containing supporting scripts. This is normally the ``scripts/``
   subdirectory in the model root directory ``modeldir``.

You are now ready to build and execute a default run. To do this, simply execute the job script
as follows::

  $ tiegcm-linux.job &

The compilation output will be displayed. If the build is successful (and exec=TRUE),
the model will be executed, and stdout will go to the specified ``output`` file.
If the job is successful, you can edit and rename the namelist input file, reset
``input`` in the job script, and re-execute the job script. If there has been
no change to the source code, it will not need to recompile, and will use the pre-existing 
executable.

.. _ibm-systems:

.. index:: ibm, aix

Running the model on IBM/AIX Platforms
--------------------------------------

.. note::

   This section contains some information that is specific to user's
   of the NCAR IBM system ``"bluefire"``. User's of other IBM systems
   may need to make adjustments for their particular environment.
   For more information about the NCAR bluefire system, see
   http://www2.cisl.ucar.edu/docs/bluefire-user-guide

The model can be built and executed on IBM platforms running AIX with
the xlf90 (mpxlf_r) compiler. You can use the same procedure described 
in the previous section, except that you use the IBM job script 
:download:`tiegcm-ibm.job <_static/tiegcm-ibm.job>` instead of the
Linux job script tiegcm-linux.job.

The IBM job script has the same user-settable shell variables as the
Linux job script, but the default settings are slightly different::

  set modeldir = tiegcm_trunk
  set execdir  = tiegcm_trunk-aix
  #set input    = tiegcm.inp
  set output   = tiegcm.out
  set make     = Make.bluefire
  set mpi      = TRUE
  set modelres = 5.0
  set debug    = FALSE
  set exec     = TRUE
  set utildir  = $modeldir/scripts

Note the ``execdir`` name, and the ``make`` file 
:download:`Make.bluefire <_static/Make.bluefire>`

Also note the special "``#BSUB``" directives at the top of the IBM
job script (descriptions in the right-hand column are for this document
only, and are not in the script itself)::

  #BSUB -J tiegcm_trunk               # Job name
  #BSUB -P 24100004                   # NCAR project number
  ##BSUB -q regular                   # regular queue (commented here)
  ##BSUB -n 32                        # number of processors (commented here)
  #BSUB -q debug                      # debug queue
  #BSUB -n 8                          # number of processors (MPI tasks)
  #BSUB -o tiegcm_trunk.%J.out        # stdout file
  #BSUB -e tiegcm_trunk.%J.out        # stderr file
  #BSUB -N
  #BSUB -u $LOGNAME@ucar.edu          # email notification address
  #BSUB -W 1:00                       # wallclock limit (6-hr max at NCAR)

These are resource settings for the Load Sharing Facility (LSF),
the batch queuing system sold by Platform Computing. The LSF is
used for scheduling jobs on the ``bluefire`` IBM system at NCAR.
This job will be submitted to the debug queue, requesting
8 processors, with a wallclock limit of 1 hour. Note the double
pound-sign "##" indicates a commented field.

To submit the IBM job script to the LSF batch system, type::

  $ bsub < tiegcm-ibm.job

Watch the progress of your LSF job with the command::

  $ bjobs

You can kill a LSF job with this command::

  $ bkill job_ID

Where ``job_ID`` is the job identifier given in the ``bjobs`` command.

For more information about the LSF, see the Wikipedia site:

http://en.wikipedia.org/wiki/Platform_LSF

or the Platform Computing site:

http://www.platform.com/workload-management/high-performance-computing/lp
