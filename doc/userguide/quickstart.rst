
.. _quickstart:

====================
QuickStart Procedure
====================

This chapter is intended to provide a tutorial-like procedure for downloading 
source code and start-up data files, building the model, and executing a short 
default run on a 64-bit Linux system, or a supercomputer (linux cluster) like
the NCAR |ncarsuper| machine.

.. index:: download

.. _download:

Downloading the model source code and required data files
---------------------------------------------------------

The model source code and related input data files may be downloaded from 
the |tgcm_version| download page (you will need to provide an email address, 
but login and password are NOT required).

:TIEGCM download page: http://www.hao.ucar.edu/modeling/tgcm/download.php

The following tar files are available:

 * Instructions (readme.download)
 * Source Code tar file (30 MB) (tiegcm2.0.tar)
 * Startup and data tar file for 5.0-deg model (tiegcm_res5.0_data.tar) (1.5 GB)
 * Startup and data tar file for 2.5-deg model (tiegcm_res2.5_data.tar) (4 GB)
  
Download the source code tar file and the data tar file for the 5.0-deg model 
to a large scratch disk on either your Linux desktop, or the NCAR supercomputer
|ncarsuper|. To download all files and make default runs at both resolutions
you will need at least 5.5 GB of disk space. Extracting these tarballs will result 
in directories with the same names. When extracting the source code, you will 
also get these default namelist input files and job scripts (or you can download
them here):

 * :download:`tiegcm-linux.job <../../scripts/tiegcm-linux.job>`: Default csh job script for 64-bit Linux desktop
 * :download:`tiegcm-ys.job <../../scripts/tiegcm-ys.job>`: Default csh job script for the NCAR supercomputer |ncarsuper|
 * :download:`tiegcm_res5.0_default.inp <../../scripts/tiegcm_res5.0_default.inp>`: Default namelist input for 5.0-degree resolution model
 * :download:`tiegcm_res2.5_default.inp <../../scripts/tiegcm_res2.5_default.inp>`: Default namelist input for 2.5-degree resolution model

Making a Default 5-deg model run
--------------------------------

The job scripts are set up to make a short (1-day) 5-degree model run (March Equinox
Solar Minimum conditions).  At this point, you should be able to simply type the 
job script name appropriate for the current machine/system on the command line.  
The job script will create an execution directory (tiegcm.exec), and build and 
execute the model there.  If successful, the stdout log will be tiegcm_res5.0.out,
and model output netCDF history files will be in the execution directory.

.. note::
  A warning for user's of previous revisions of TIEGCM: Do not use old namelist
  input files or job scripts from previous revisions. Copy the default files
  from the :term:`scripts/` directory, and modify them for your own runs.
  Also, for initial runs, do not build/execute the model in an old :term:`execdir`. 
  Instead, allow the job script to make a new execdir for you.

Switching to 2.5-degree Model Resolution
----------------------------------------

To make a default run of the 2.5-deg model, edit the job script and reset 4 shell
variables as follows:

 * set tgcmdata = tiegcm_res2.5_data
 * set input    = tiegcm_res2.5.inp
 * set output   = tiegcm_res2.5.out
 * set modelres = 2.5

If you are on the NCAR supercomputer |ncarsuper|, you should also make the
following changes to tiegcm-ys.job, to use 64 cores:

 * #BSUB -n 64
 * #BSUB -R "span[ptile=16]"

Then execute the job script for the default 2.5-deg model run.

.. _continuation_run:

Making a Continuation (Restart) Run
-----------------------------------

A model run can be continued (restarted) from the last primary history written
by the previous run. To do this, you must modify the namelist input file
as follows (refer to :ref:`namelist_params` for more information):

  1. Comment or remove SOURCE and SOURCE_START if the previous run was an initial run.
  2. Reset START_DAY, START and STOP as necessary
  3. Make sure one of the files in the OUTPUT list contains the new START history
  4. Increment the starting volume number of SECOUT (pre-existing secondary
     output history files will be overwritten).

Moving to "Production" Mode
---------------------------

When you are ready to make longer model runs, or especially if you
are planning to modify the source code, its best to move your working
directory (with the model source directory, and any job scripts or
namelist input files) to a disk space that is regularly backed up, 
e.g., under your home. You can leave the data (:term:`tgcmdata`)
and execution directories (:term:`execdir`) on the large scratch disk, 
but you must then set the tgcmdata and execdir to absolute paths in
the job script. 

As you proceed, you can create new working directories, and corresponding 
execdirs as needed. If you modify the source code, the job script will 
call gmake, and dependent source files will be recompiled as necessary. 
If you switch between resolutions using the same execdir, the entire code 
will be recompiled for the new resolution.

Notes for Users on the HAO network
----------------------------------

 * Startup and data files for |tgcm_version| are in /hao/aim/tgcm/data/tiegcm2.0
 * The /hao/aim disk can be slow so its best (and probably fastest) to run the 
   model on the local Linux desktop disk (e.g., something like: set 
   execdir = /export/data1/$user/tiegcm.exec)
 * Its usually best to run with 4 or 8 processors on the Linux box (set nproc = 4).
 * Although the model has been built with PGI and gfortran at hao, the model will
   run fastest if built with the Intel compiler (set make = Make.intel_hao64)

Notes for Users on the NCAR /glade disk (|ncarsuper|)
-----------------------------------------------------

 * For more information on using |ncarsuper|, see 
   `NCAR CISL documentation <http://www2.cisl.ucar.edu/resources/computational-systems/yellowstone>`_
 * Startup and data files for |tgcm_version| are in /glade/p/hao/tgcm/data/tiegcm2.0

The |ncarsuper| system uses the Load Sharing Facility (LSF) as a batch job management system:

 * See `LSF Introduction Guide <http://www.vub.ac.be/BFUCC/LSF/>`_ for a brief overview.
 * Also see `CISL Platform LSF job script examples <https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/using-computing-resources/running-jobs/platform-lsf-job-script-examples>`_
 * LSF resources are specified using LSF #BSUB commands. The default |ncarsuper| job script tiegcm-ys.job
   uses the following settings::

   #BSUB -J tiegcm                 [arbitrary job name]
   #BSUB -P P28100036              [your authorized NCAR project number (this one is used at hao)]
   #BSUB -q premium                [queue (can be regular, premium, standby, etc)]
   #BSUB -o tiegcm.%J.out          [specify stdout file (different from model stdout)]
   #BSUB -e tiegcm.%J.out          [specify stderr]
   #BSUB -N                        [not sure what this is for]
   #BSUB -u $LOGNAME@ucar.edu      [send email to this address after job has completed]
   #BSUB -W 0:30                   [wallclock limit hh:mm (max 12 hours at NCAR)]
   #BSUB -n 16                     [use 16 processors (64 for 2.5-deg tiegcm)]

 * To calculate wallclock time for a 5-deg run with 16 cores, use .07 secs/timestep.
   For example, a 1-day simulation with a 60 sec timestep: ((24*3600) / 60 * .07) / 60 = 1.68 minutes

 * To calculate wallclock time for a 2.5-deg run with 64 cores, use .15 secs/timestep. 
   For example, a 1-day simulation with a 30 sec timestep: ((24*3600) / 30 * .15) / 60 = 7.2 minutes
     
