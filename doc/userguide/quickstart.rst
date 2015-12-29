
.. _quickstart:

====================
QuickStart Procedure
====================

.. note::

   This document is up to date for version |version| of the TIEGCM

This document is intended to provide a tutorial-like procedure for downloading 
source code and start-up data files, building the model, and executing a short 
default run on a 64-bit Linux system.

.. index:: download

.. _download:

Downloading the model source code and required data files
---------------------------------------------------------

The model source code and related input data files may be downloaded from 
the TIEGCM download page of the main TGCM website:

:download_url:`download.php`

You will need to provide an email address (login and password are NOT required).
Documentation and Postprocessor codes are also available on the download site,
but all you need for now is the source code, and corresponding data files.
Both are provided as gzipped tar files. There are separate data files for 
5.0-degree and 2.5-degree model resolutions.

Download the model source code:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Make an empty working directory to contain the model source code (:term:`/home/user/tiegcm_wrkdir`)
#. Download the file |tgcm_version|.tar.gz from the :download_url:`download <download.php>` page to the working directory.
#. Uncompress (gunzip) and extract the tar file.

At this point, you should have something like this in your working directory::

  total 64
  -rw-r--r--  1 foster  staff  5162 Dec 15 15:39 README.download
  -rw-r--r--  1 foster  staff  5608 Dec 15 15:39 Release_Notes
  -rwxr-xr-x  1 foster  staff  6713 Dec 15 15:39 tiegcm-linux.job*
  drwxr-xr-x  8 foster  staff   272 Dec 15 13:21 tiegcmx.x/
  -rw-r--r--  1 foster  staff  6116 Dec 15 15:39 tiegcmlicense.txt
 
These files and directories contain the following:

.. describe:: README.download

   :download:`README.download <_static/README.download>`
   Instructions for building and making a short default run.

.. describe:: Release_Notes

   Release notes for this version of the model.

.. describe:: tiegcm-linux.job 

   Job script for building and executing under Linux (64-bit) systems.
   (:download:`tiegcm-linux.job <_static/tiegcm-linux.job>`)

.. describe:: tiegcmx.x/

   Model root directory for this version (e.g., |tgcm_version|), containing source code, supporting scripts, and documentation.

.. describe:: tiegcmlicense.txt 

   Academic license agreement.
   (:download:`tiegcmlicense.txt <_static/tiegcmlicense.txt>`)

Download required data files:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TGCMDATA

#. Make a directory on a large temporary storage disk (e.g., :term:`/mytmpdir/tiegcm/data`).
#. Set environment variable :term:`TGCMDATA` to the absolute path of this directory.
#. Download the file |download_5.0data| from the :download_url:`download <download.php>` 
   page to the data directory.
#. Uncompress (gunzip) and extract the tar file.

Now, you should have something like this in your data directory (:term:`TGCMDATA`)::

  -rw-r--r-- 1 foster hao  1695912 Nov  4 14:18 gpi_1960001-2015090.nc
  -rw-r--r-- 1 foster hao 23889780 Dec  7 09:39 gswm_diurn_5.0d_99km.nc
  -rw-r--r-- 1 foster hao 23889816 Dec  7 09:39 gswm_nonmig_diurn_5.0d_99km.nc
  -rw-r--r-- 1 foster hao 23889768 Dec  7 09:39 gswm_nonmig_semi_5.0d_99km.nc
  -rw-r--r-- 1 foster hao 23889784 Dec  7 09:39 gswm_semi_5.0d_99km.nc
  -rw-r--r-- 1 foster hao  5573432 Nov  4 14:09 imf_OMNI_2003001-2003365.nc
  -rw-r--r-- 1 foster hao  5573432 Nov  4 14:56 imf_OMNI_2006001-2006365.nc
  -rw-r--r-- 1 foster hao  5588648 Nov  4 21:45 imf_OMNI_2008001-2008366.nc
  -rw-r--r-- 1 foster hao 16990860 Nov 21 20:06 tiegcm_res5.0_climatology_smax_prim.nc
  -rw-r--r-- 1 foster hao 16990860 Nov 21 20:06 tiegcm_res5.0_climatology_smin_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_dec2006_heelis_gpi_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_dec2006_weimer_imf_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_decsol_smax_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_decsol_smin_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_junsol_smax_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_junsol_smin_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_mareqx_smax_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_mareqx_smin_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_nov2003_heelis_gpi_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_nov2003_weimer_imf_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_sepeqx_smax_prim.nc
  -rw-r--r-- 1 foster hao 16990848 Nov 12 14:04 tiegcm_res5.0_sepeqx_smin_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_whi2008_heelis_gpi_prim.nc
  -rw-r--r-- 1 foster hao 16990868 Nov 12 14:04 tiegcm_res5.0_whi2008_weimer_imf_prim.nc

These are all netCDF data files (if you have netcdf installed on your system, you can
see the contents of these files with the command "ncdump -c file.nc").

.. describe:: gpi*.nc

  | Files containing dated geophysical indices and solar fluxes Kp, f107, and f107a 
  | Specify these files with the namelist input keyword :ref:`GPI_NCFILE <GPI_NCFILE>`

.. describe:: gswm*.nc (GSWM* namelist input keywords)

  | Lower boundary tidal perturbations from the Global Scale Wave Model `GSWM <http://www.hao.ucar.edu/modeling/gswm/gswm.html>`_
  | Specify these files with the namelist input keywords :ref:`GSWM data files <GSWM>`

.. describe:: imf*.nc 

  | Interplanetary Magnetic Field OMNI data files
  | For use when namelist input keyword :ref:`POTENTIAL_MODEL <POTENTIAL_MODEL>` = 'WEIMER'
  | Specify these files with the namelist input keyword :ref:`IMF_NCFILE <IMF_NCFILE>`

.. describe:: tiegcm_res5.0*.nc

  | Startup history files for an initial run (of tiegcm in this case).
  | Specify these files with the namelist input keyword :ref:`SOURCE <SOURCE>`

.. index:: resolution
| For 5.0-degree model :term:`resolution`, the data files total about 500 MB.
| For 2.5-degree model :term:`resolution`, the data files total about 2.6 GB.

Making a default run on a 64-bit Linux system
---------------------------------------------

.. _jobscript:


Take a look at the Linux job script 
:download:`tiegcm-linux.job <_static/tiegcm-linux.job>`. Near the top are
several shell variables, with their default settings, which configure the 
job script (variables and values may vary between model versions)::

  set modeldir = tiegcmx.x
  set execdir  = /hao/aim/$user/tiegcm-linux
  set input    = $modeldir/scripts/tiegcm_res5.0_default.inp
  set output   = tiegcm.out
  set make     = Make.intel_hao64
  set modelres = 5.0
  set mpi      = TRUE
  set nproc    = 4
  set debug    = FALSE
  set exec     = TRUE
  set utildir  = $modeldir/scripts

Following are brief explanations of the job script shell variables:

.. note::
   
   Absolute or relative paths are acceptable when specifying directories.
   Relative paths should be relative to the *working directory* (:term:`workdir`).

.. describe:: modeldir

   The model root directory (:term:`modeldir` from the source code download). 
   This will contain subdirectories :term:`src/` , :term:`scripts/` , :term:`doc/` , 
   and :term:`tgcmrun/` .

.. describe:: execdir

   This is the execution directory (:term:`execdir`), in which the model will
   be built and executed. It will be created if it does not already exist.
   This directory will also contain the model output netCDF history files.

.. describe:: input

   The :ref:`namelist input file <namelist>`. The default namelist file is in the scripts directory
   under the model root with file name :download:`tiegcm_res5.0_default.inp <_static/tiegcm_res5.0_default.inp>`
   (for 5-degree resolution), or :download:`tiegcm_res2.5_default.inp <_static/tiegcm_res2.5_default.inp>`
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
    * :download:`Make.intel_hao64 <_static/Make.intel_hao64>` (for Intel compiler)
    * :download:`Make.pgi_hao64 <_static/Make.pgi_hao64>` (for PGI compiler)
    * :download:`Make.gfort_hao64 <_static/Make.gfort_hao64>` (for gfortran compiler)

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

Running the model on NCAR Linux Supercomputer |ncarsuper|
---------------------------------------------------------

.. note::

   This section contains information that is specific to user's
   of the NCAR Linux Supercomputer |ncarsuper|:

   | Linux yslogin3 2.6.32-358.el6.x86_64 #1 SMP 
   | Tue Jan 29 11:47:41 EST 2013 x86_64 x86_64 x86_64 GNU/Linux

   For more information about the NCAR |ncarsuper| system, see
   |ncarsuper_url|

The model can be built and executed on |ncarsuper|, using the Intel
compiler and the intelmpi implementation. To do this, copy and modify
the job script tiegcm-ys.job from the scripts directory.

The |ncarsuper| :download:`tiegcm-ys.job <_static/tiegcm-ys.job>` has 
the same user-settable shell variables as the Linux job script, but 
the default settings are slightly different::

  set modeldir  = tiegcmx.x
  set execdir   = /glade/scratch/$user/tiegcm_trunk/tiegcm-ys
  set input     = $modeldir/scripts/tiegcm_res5.0_default.inp
  set output    = tiegcm.out
  set make      = Make.intel_ys
  set mpi       = TRUE
  set modelres  = 5.0
  set debug     = FALSE
  set exec      = TRUE
  set utildir   = $modeldir/scripts
  set runscript = run.lsf

Where ``x.x`` refers to the version number (|version|).
Note the :term:`execdir` name, and the ``make`` file 
:download:`Make.intel_ys <_static/Make.intel_ys>`.
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
