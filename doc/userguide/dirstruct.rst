
.. _dirstruct:

Directory Structure
===================

To make a model run, some directory paths must be defined:

* :term:`workdir`: User created directory from which model runs are submitted.

* :term:`modeldir`: Directory containing the model source code (may be an svn working copy)

* :term:`tgcmdata`: Location of start-up history and input data files (see :ref:`datadir` below)

* :term:`execdir`: Directory in which the model will be built and executed 

These paths are set by the user in the job script, as in these examples:
(these can be absolute paths or relative to the working directory)::

  set modeldir = tiegcm_trunk
  set tgcmdata = /my/prettybig/data/tiegcm.data
  set execdir  = /my/big/exec/disk/tiegcm.exec

See this :download:`example job script <../../scripts/tiegcm-linux.job>`

.. _workdir:

Working Directory (created by the user)
---------------------------------------

To get started, the user will typically create a working directory from
which model runs are submitted. The user's working directory will typically 
look something like this::

                      workdir
                         |
 -----------------------------------------------
              |                    |
            *.inp              modeldir/             
            *.job                          
            *.out                         
                                
Here, \*.inp are :term:`namelist input` files, \*.job are 
:term:`job script`'s, and \*.out are stdout :term:`output` files from model 
and are validated by the input module (input.F). 

The job script in your working directory contains a shell variable specifying
the path to the :term:`modeldir`, so it knows where to find the source code and 
supporting scripts for the build process. The namelist input file also refers to 
the :term:`datadir` path for start-up and other data input files (e.g., :ref:`SOURCE <SOURCE>`, 
:ref:`GPI_NCFILE <GPI_NCFILE>`, :ref:`IMF_NCFILE <IMF_NCFILE>`, etc). 
These namelist parameters can use the environment variable :term:`TGCMDATA` to 
specify the :term:`datadir` (see section on :ref:`namelist input files <namelist>`).
The job script shell variable :term:`tgcmdata`, if set, will override the TGCMDATA env var.

.. _modeldir:

Model Directory (model source code and supporting utilities)
------------------------------------------------------------

The model root directory is what you get when you :ref:`download <download>` the 
model source code tar file, or check out the code from the svn repository. 
The model directory contains subdirectories with the model source code, 
supporting scripts, documentation, and a python code to make test and benchmark runs::

                                modeldir
                                   |
   -------------------------------------------------------------------------
      |               |                |                 |             |
     src/          scripts/           doc/            tgcmrun/     benchmarks/
      |               |                |                 |             |
     *.F90          Make.*         userguide/           *.py      run_climatology
     *.F          linux.job       description/          run_*       run_seasons
     *.h           ibm.job          release/          tgcmrun       run_storms
                 default.inp       diags.table                    archive_hpss 
                tgcm_contents      perf.table                     make_listings
                 tgcm_ncdump    README.download                     postproc/
                    etc                           

:term:`src/` directory contents:

* Fortran source code \*.F, \*.F90, \*.h. The source code is f90 standard compliant, and most 
  source files are in fixed-format fortran. There is a single header file, defs.h,
  which contains grid definitions and dimensions for different :term:`resolution` s.

:term:`scripts/` directory contents:

* **Make.\***: Makefiles containing platform-dependent compiler flags, 
  Make variables, and library locations.  These files be copied, renamed, 
  and customized for the user's platform/machine environment. 
* **Make.*_hao64**: Three compilers are supported on the linux desktop platform: 
  intel, pgi, and gfortran.
* **Make.intel_ys**:  Makefile for intel compiler on the NCAR supercomputer |ncarsuper|.
* **Makefile**: The main makefile. The Make.xxx file currently in use is included in the Makefile 
  at build time.
* **tiegcm-linux.job**: Default model :term:`job script` for Linux desktop systems.
* **tiegcm-ys.job**: Default model :term:`job script` for the NCAR |ncarsuper| supercomputer.
* **tiegcm_res5.0_default.inp**: Default namelist input file for 5.0-degree resolution.
* **tiegcm_res2.5_default.inp**: Default namelist input file for 2.5-degree resolution.
* **download**: Directory in which to make source and data tar files for :ref:`download <download>` from the TGCM website

There are several additional utilities in the scripts directory that are used by
the build system or by the user to perform various tasks or to obtain information
(see :download:`README in scripts directory <../../scripts/README>` for more information).
directory for more information.

:term:`doc/` directory contents:

* **userguide/**: Directory containing `Python Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_ source files for the User's Guide (this document)

* **description/**: Directory containing source files for the 
  `Model Description <http://www.hao.ucar.edu/modeling/tgcm/doc/description/model_description.pdf>`_

* **release/**: Directory containing source files for the 
  :base_url:`Release Documentation <release/html>`

* **diags.table**: :download:`Table of diagnostic fields <_static/diags.table>` that can be 
  saved on secondary history files.

.. index:: perf.table

* **perf.table**: :download:`Table of performance statistics <_static/perf.table>` for both
  models (tiegcm and timegcm) at both :term:`resolution`.

* **README.download**: :download:`Instructions <../../scripts/download/README.download>` for how to make a 
  quick-start default build and execution of the model after downloading the source code and data.

:term:`tgcmrun/` directory contents:

* Python code to make benchmark runs interactively or from shell scripts. Type 'tgcmrun' on 
  the command line for interactive, or execute the run_* scripts to make benchmark series runs.
* For more information on benchmark runs made for the current release, please see 
  :base_url:`Release Documentation <release/html>`

:term:`benchmarks/` directory contents:

* Shell scripts that call :term:`tgcmrun/` to make benchmark runs:

  * run_climatology
  * run_seasons
  * run_storms
  * run_perf

* Script archive_hpss to archive benchmark runs on the hpss 
  (see `HPSS <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/hpss>`_)
* Script make_listings for making lists of files related to benchmark runs
* Subdirectory postproc/ contains scripts that call the tgcmproc utility to post-process
  benchmark runs.

.. _datadir:

Data Directory (startup and data files)
---------------------------------------

The public |modeluc| data directory is what you get when you :ref:`download <download>` 
the data tar file. This directory is typically referred to with the environment variable
:term:`TGCMDATA`, but can be set with the :term:`tgcmdata` shell variable in the job script
(the shell variable, if set, will take precedence over the environment variable)::

                  datadir for tiegcmx.x
                          |
   ----------------------------------------------
                          |
                tiegcmx.x_res5.0_*.nc
                tiegcmx.x_res2.5_*.nc 
                        gpi*.nc
                      gswm*5.0d*.nc
                      gswm*2.5d*.nc
                      imf_OMNI_*.nc
                         etc

These are netCDF history startup and data files for running the current version of the
model ( |tgcm_version| )
They are specified in the namelist input file (see :ref:`namelist input files <namelist>` 
for more information). These files are available for download, see :ref:`download`.


* **tiegcmx.x_res5.0_*.nc**: History start-up files for the 5.0-degree resolution
  model. These files contain a single history with initial conditions for starting 
  the model at a specific date and time. These are typically the first history
  for a benchmark run (seasonal, storm simulations, and climatologies). 
  Namelist input parameter: :ref:`SOURCE <SOURCE>`. Here are the files for |tgcm_version|:

  | tiegcm_res5.0_climatology_smax_prim.nc
  | tiegcm_res5.0_climatology_smin_prim.nc
  | tiegcm_res5.0_dec2006_heelis_gpi_prim.nc
  | tiegcm_res5.0_dec2006_weimer_imf_prim.nc
  | tiegcm_res5.0_decsol_smax_prim.nc
  | tiegcm_res5.0_decsol_smin_prim.nc
  | tiegcm_res5.0_junsol_smax_prim.nc
  | tiegcm_res5.0_junsol_smin_prim.nc
  | tiegcm_res5.0_mareqx_smax_prim.nc
  | tiegcm_res5.0_mareqx_smin_prim.nc
  | tiegcm_res5.0_nov2003_heelis_gpi_prim.nc
  | tiegcm_res5.0_nov2003_weimer_imf_prim.nc
  | tiegcm_res5.0_sepeqx_smax_prim.nc
  | tiegcm_res5.0_sepeqx_smin_prim.nc
  | tiegcm_res5.0_whi2008_heelis_gpi_prim.nc
  | tiegcm_res5.0_whi2008_weimer_imf_prim.nc

* **tiegcmx.x_res2.5_*.nc**: History start-up files for the 2.5-degree resolution
  model. These files contain a single history with initial conditions for starting 
  the model at a specific model date and time. These are typically the first history
  for a benchmark run (seasonal, storm simulations, and climatologies).
  Namelist input parameter: :ref:`SOURCE <SOURCE>` Here are the files for |tgcm_version|:

  | tiegcm_res2.5_bgndlbc_hwm_msis.nc
  | tiegcm_res2.5_bgndlbc_saber_hrdi.nc
  | tiegcm_res2.5_climatology_smax_prim.nc
  | tiegcm_res2.5_climatology_smin_prim.nc
  | tiegcm_res2.5_dec2006_heelis_gpi_prim.nc
  | tiegcm_res2.5_dec2006_weimer_imf_prim.nc
  | tiegcm_res2.5_decsol_smax_prim.nc
  | tiegcm_res2.5_decsol_smin_prim.nc
  | tiegcm_res2.5_junsol_smax_prim.nc
  | tiegcm_res2.5_junsol_smin_prim.nc
  | tiegcm_res2.5_mareqx_smax_prim.nc
  | tiegcm_res2.5_mareqx_smin_prim.nc
  | tiegcm_res2.5_nov2003_heelis_gpi_prim.nc
  | tiegcm_res2.5_nov2003_weimer_imf_prim.nc
  | tiegcm_res2.5_sepeqx_smax_prim.nc
  | tiegcm_res2.5_sepeqx_smin_prim.nc
  | tiegcm_res2.5_whi2008_heelis_gpi_prim.nc
  | tiegcm_res2.5_whi2008_weimer_imf_prim.nc
 
* **gpi\*.nc**
  GeoPhysical Indices data files (3-hourly Kp and F10.7 cm solar flux).
  Namelist Input parameter: :ref:`GPI_NCFILE <GPI_NCFILE>`
  
* **gswm*5.0d*.nc**
  `Global Scale Wave Model <http://www.hao.ucar.edu/modeling/gswm/gswm.html>`_
  data files, used to specify tidal perturbations for the lower boundary of
  the |modeluc| for the 5-degree resolution. There are 4 separate files for 
  diurnal, semi-diurnal, migrating and non-migrating tides. 
  Namelist Input parameter: :ref:`GSWM <GSWM>`.
  
* **gswm*2.5d*.nc**
  `Global Scale Wave Model <http://www.hao.ucar.edu/modeling/gswm/gswm.html>`_
  data files, used to specify tidal perturbations for the lower boundary of
  the |modeluc| for the 2.5-degree resolution. There are 4 separate files for 
  diurnal, semi-diurnal, migrating and non-migrating tides. 
  Namelist Input parameter: :ref:`GSWM <GSWM>`.

* **imf_OMNI_*.nc**
  Interplanetary Magnetic Field OMNI data files. Namelist read parameter is
  :ref:`IMF_NCFILE <IMF_NCFILE>`. These files contain data for the BX,BY,BZ 
  components of the IMF, solar wind velocity and solar wind density.
  See `HAO public ftp page <http://download.hao.ucar.edu/pub/tgcm/data>`_ 
  to download imf data files for years not included on the |tgcm_version|
  data download.

.. _execdir:

Execution Directory (build and execute)
---------------------------------------

The model is built and executed in the execution directory (:term:`execdir`). 
The path to the execution directory is specified by the execdir shell variable
in the :term:`job script`. The job script will create the execdir for you if
it does not already exist. The following file types are typically found in the 
execution directory:

.. note::
  When making your first run, its best to let the job script create the execdir
  for you. It is not wise to use an execdir used for revisions prior to |tgcm_version|.
  Also, if you have build or execution problems, it will sometimes help to remove
  the execdir and let the job script start over.

* ***.o**:
  Object files produced by the compiler.

* ***.mod**:
  Module files produced by the compiler.

* ***PET*LogFile**:
  ESMF log files.

* **tiegcm*.nc**:
  Model output netCDF history files.

* **M***:
  Makefiles.

The model executable also resides in the execution directory.
