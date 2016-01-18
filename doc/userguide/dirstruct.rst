
.. _dirstruct:

Directory Structure
===================

This chapter describes a typical directory structure the user will be
working with when running the |modeluc| model. The working directory
(:term:`workdir`) is the "root" directory of the user's project. 

The model directory (:term:`modeldir`) is typically a subdirectory under
the working directory, and contains the model source code, supporting
scripts, documentation, and a semi-automated validation and testing program.

The data directory (:term:`datadir`) is typically on a large temporary disk 
capable of storing start-up and input data files necessary for running the model.
The execution directory (:term:`execdir`) is also on a large temporary disk,
and contains the compiled object and module files, the model executable, and 
the model output netCDF history files.

.. _workdir:

Working Directory (:term:`workdir`)
-----------------------------------

The user's working directory will typically look something like this
(the :term:`datadir` can be on a large separate disk system)::

                      workdir
                         |
 -----------------------------------------------
              |                    |
            *.inp              modeldir/             
            *.job                          
            *.out                         
                                
Here, *.inp are :term:`namelist input` files, *.job are 
:term:`job script` s, and *.out are stdout :term:`output` files from model 
runs. 

The job script(s) in your working directory contains a shell variable specifying
the path to the :term:`modeldir`, so it knows where to find the source code and 
supporting scripts for the build process. The namelist input file also refers to 
the :term:`datadir` path for start-up and other data input files (e.g., :ref:`SOURCE <SOURCE>`, 
:ref:`GPI_NCFILE <GPI_NCFILE>`, :ref:`IMF_NCFILE <IMF_NCFILE>`, etc). 
These namelist parameters can use the environment variable :term:`TGCMDATA` to 
specify the :term:`datadir` (see section on :ref:`namelist input files <namelist>`).

.. _modeldir:

Model Directory (:term:`modeldir`)
----------------------------------

The model root directory is what you get when you :ref:`download <download>` the 
model source code tar file. The model directory contains subdirectories with the 
model source code, supporting scripts, documentation, and a python code to make
test and benchmark runs::

                                modeldir
                                   |
   ----------------------------------------------------------------------
      |               |                |                 |
     src/          scripts/           doc/            tgcmrun/
      |               |                |                 | 
     *.F90          Make.*         userguide/           *.py 
     *.F          linux.job       description/          run_* 
     *.h           ibm.job          release/          tgcmrun 
                 default.inp       diags.table       
                tgcm_contents      perf.table
                 tgcm_ncdump    README.download    
                    etc                           

:term:`src/` directory contents:

* Fortran source code *.F, *.F90, *.h. The source code is f90 standard compliant, and most 
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

There are several additional utilities in the scripts directory that are used by
the build system or by the user to perform various tasks or to obtain information
(see :download:`README.scripts <_static/README.scripts>` for more information).
directory for more information.

:term:`doc/` directory contents:

* **userguide/**: Directory containing `Python Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_ source files for the User's Guide (this document)

* **description/**: Directory containing source files for the 
  :base_url:`TIEGCM Model Description <description/tiegcm_modeldes_6Oct09.pdf>`

* **release/**: Directory containing source files for the 
  :base_url:`Release Documentation <release/html>`

* **diags.table**: :download:`Table of diagnostic fields <_static/diags.table>` that can be 
  saved on secondary history files.

.. index:: perf.table

* **perf.table**: :download:`Table of performance statistics <_static/perf.table>` for both
  models (tiegcm and timegcm) at both :term:`resolution`.

* **README.download**: :download:`Instructions <_static/README.download>` for how to make a 
  quick-start default build and execution of the model after downloading the source code and data.

:term:`tgcmrun/` directory contents:

* Python code to make benchmark runs interactively or from shell scripts. Type 'tgcmrun' on 
  the command line for interactive, or execute the run_* scripts to make benchmark series runs.
* For more information on benchmark runs made for the current release, please see 
  :base_url:`Release Documentation <release/html>`

.. _datadir:

Data Directory (:term:`datadir`)
--------------------------------

The public |modeluc| data directory is what you get when you :ref:`download <download>` 
the data tar file. This directory is typically referred to with the environment variable
:term:`TGCMDATA`.:: 

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

These are netCDF history startup and data files for running the current version of the
model ( |tgcm_version| )
They are specified in the namelist input file (see :ref:`namelist input files <namelist>` 
for more information). Additional files may be downloaded from the 
:ftp_url:`HAO public ftp page <>`


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
  See :ftp_url:`HAO public ftp page <>` to download imf data files for years 
  not included at the :download_url:`tgcm download website <>`.

.. _execdir:

Execution Directory (:term:`execdir`)
-------------------------------------

The model is built and executed in the execution directory (:term:`execdir`). 
The path to the execution directory is specified by the execdir shell variable
in the :ref:`job script <jobscript>`. The following file types are typically
found in the execution directory:

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
